"""Publish a dated, non-vacuous snapshot of completed numerical controls."""
import run_transient as rt
import json,shutil
from datetime import datetime,timezone,timedelta
import numpy as np


def main(root=None):
    root=rt.ROOT/'output/transient_converged' if root is None else rt.Path(root);rt.configure_scratch(root)
    rows=[];band_rows=[]
    for kind,mode,index in [('guided','A0',1),('bulk','P',0)]:
        folder=root/'spatial_checks'/kind
        profile=json.loads((folder/'base/group/settings.json').read_text())['layout']['profile']
        base=np.diag(np.load(folder/'base/group/spectra.npz')['T'][0])/np.diag(np.load(folder/'base/reference/spectra.npz')['T'][0])
        for variant in (['mesh','absorber'] if kind=='guided' else ['mesh','absorber','width']):
            row=dict(kind=kind,mode=mode,profile=profile,check=variant,status='pending',center_amplitude_drop_percent=float(100*(1-abs(base[index]))))
            if all((folder/variant/c/'spectra.npz').exists() for c in ['reference','group']):
                other=np.diag(np.load(folder/variant/'group/spectra.npz')['T'][0])/np.diag(np.load(folder/variant/'reference/spectra.npz')['T'][0])
                delta=float(abs(abs(base[index])-abs(other[index])))
                phase=float(abs(np.angle(base[index]/other[index]))*180/np.pi)
                relative=delta/max(abs(1-abs(base[index])),abs(1-abs(other[index])),1e-12)
                row.update(magnitude_delta_percentage_points=delta*100,phase_error_degrees=phase,relative_loss_change_percent=relative*100,
                    status='passed' if delta<=2e-4 and phase<=1 and relative<=.05 else 'failed')
            rows.append(row)
            band_row=dict(kind=kind,mode=mode,profile=profile,check=variant,status='pending')
            if all((folder/v/c/'spectra_band.npz').exists() for v in ['base',variant] for c in ['reference','group']):
                def bands(v,c):return np.load(folder/v/c/'spectra_band.npz')['T'][:,index,index]
                aa=bands('base','group')/bands('base','reference');bb=bands(variant,'group')/bands(variant,'reference')
                delta=float(np.max(abs(abs(aa)-abs(bb))));phase=float(np.max(abs(np.angle(aa/bb)))*180/np.pi)
                relative=float(np.max(abs(abs(aa)-abs(bb))/np.maximum(np.maximum(abs(1-abs(aa)),abs(1-abs(bb))),1e-12)))
                band_row.update(magnitude_delta_percentage_points=delta*100,phase_error_degrees=phase,relative_loss_change_percent=relative*100,
                    status='passed' if delta<=2e-4 and phase<=1 and relative<=.05 else 'failed')
            band_rows.append(band_row)
    gate_path=root/'validation.json';gate=json.loads(gate_path.read_text()) if gate_path.exists() else {}
    guided_controls=[c for c in gate.get('checks',[]) if c['kind']=='guided' and c['case'] in ['reference','group']]
    guided_passed=len(guided_controls)==2 and all(c['passed'] for c in guided_controls)
    stamp=datetime.now(timezone(timedelta(hours=8))).isoformat(timespec='seconds')
    payload=dict(updated_at=stamp,scope='Completed spatial controls and separate transient acceptance; missing checks remain pending.',
        ready_for_ensemble=bool(gate.get('ready_for_ensemble',False)),guided_transient_passed=guided_passed,transient_checks=gate.get('checks',[]),checks=rows,frequency_band_checks=band_rows,C_free_gib=shutil.disk_usage('C:/').free/2**30,
        E_free_gib=shutil.disk_usage('E:/').free/2**30)
    rt.write_json(root/'LATEST_VALIDATION.json',payload)
    names={'mesh':'网格加密','absorber':'吸收区变化','width':'周期单元复制'}
    statuses={'passed':'通过','failed':'未通过','pending':'待完成'}
    table=[]
    for row in rows:
        table.append([row['mode'],names[row['check']],f'{row["magnitude_delta_percentage_points"]:.6g}' if row['status']!='pending' else '—',
                      f'{row["relative_loss_change_percent"]:.4g}%' if row['status']!='pending' else '—',
                      f'{row["phase_error_degrees"]:.4g}°' if row['status']!='pending' else '—',statuses[row['status']]])
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans'];plt.rcParams['axes.unicode_minus']=False
    fig,ax=plt.subplots(figsize=(11,3.8));ax.axis('off')
    headers=['模态','检查','幅值差（百分点）','损失相对变化','相位差','中心频率验收']
    grid=ax.table(cellText=table,colLabels=headers,cellLoc='center',loc='center',colWidths=[.08,.19,.22,.19,.14,.18])
    grid.auto_set_font_size(False);grid.set_fontsize(10);grid.scale(1,1.8)
    for i,row in enumerate(rows,1):grid[i,5].set_facecolor({'passed':'#d9eddd','failed':'#f5d8d8','pending':'#eeeeee'}[row['status']])
    fig.suptitle('最新中心频率验证：完整收敛验收尚未完成' if not payload['ready_for_ensemble'] else '数值收敛验收已通过',fontsize=14)
    fig.text(.5,.07,'门限：幅值差 ≤ 0.02 个百分点；非零损失相对变化 ≤ 5%；相位差 ≤ 1°。',ha='center',fontsize=10)
    fig.text(.5,.015,'周期复制检查不等同于随机代表体积收敛。完整验收还包含频带、时间步、观察窗和瞬态／频域一致性。',ha='center',fontsize=9)
    fig.tight_layout(rect=(.01,.13,.99,.93));fig.savefig(root/'latest_center_validation.png',dpi=160);plt.close(fig)
    band_table=[]
    for row in band_rows:
        band_table.append([row['mode'],names[row['check']],f'{row["magnitude_delta_percentage_points"]:.6g}' if row['status']!='pending' else '—',
            f'{row["relative_loss_change_percent"]:.4g}%' if row['status']!='pending' else '—',f'{row["phase_error_degrees"]:.4g}°' if row['status']!='pending' else '—',statuses[row['status']]])
    fig,ax=plt.subplots(figsize=(11,3.8));ax.axis('off')
    band_headers=headers[:-1]+['17 点频带验收']
    grid=ax.table(cellText=band_table,colLabels=band_headers,cellLoc='center',loc='center',colWidths=[.08,.19,.22,.19,.14,.18])
    grid.auto_set_font_size(False);grid.set_fontsize(10);grid.scale(1,1.8)
    for i,row in enumerate(band_rows,1):grid[i,5].set_facecolor({'passed':'#d9eddd','failed':'#f5d8d8','pending':'#eeeeee'}[row['status']])
    fig.suptitle('激励频带空间验证：各列为 17 个频点中的最大变化',fontsize=14)
    fig.text(.5,.07,'A0：131.25–268.75 kHz；P：3.28125–6.71875 MHz。',ha='center',fontsize=10)
    fig.text(.5,.015,'空间验证通过后仍需通过时间步、观察窗及瞬态／频域一致性检查，才能开展统计扩展。',ha='center',fontsize=9)
    fig.tight_layout(rect=(.01,.13,.99,.93));fig.savefig(root/'latest_band_validation.png',dpi=160);plt.close(fig)
    lines=['# 最新验证快照','',f'更新时间：{stamp}','',
        '**完整数值验收尚未全部通过；长度和随机样本扩展未启动。**' if not payload['ready_for_ensemble'] else '**完整数值验收已通过。**','',
        '|'+ '|'.join(headers)+'|','|---|---|---:|---:|---:|---|']
    lines += ['|'+'|'.join(row)+'|' for row in table]
    lines += ['','## 激励频带检查','', '|'+ '|'.join(band_headers)+'|','|---|---|---:|---:|---:|---|']
    lines += ['|'+'|'.join(row)+'|' for row in band_table]
    lines += ['', '## 瞬态验收', '', '|波型／算例|时间步减半|瞬态与频域：幅值差（百分点）|非零损失相对差|时间窗复信号差|完整验收|',
              '|---|---|---:|---:|---:|---|']
    for kind in ['guided','bulk']:
        for case in ['reference','group']:
            record=next((c for c in gate.get('checks',[]) if c['kind']==kind and c['case']==case),None)
            if record is None:
                lines.append(f'|{kind}／{case}|待完成|—|—|—|待完成|');continue
            time_control=record['differences'].get('time',{}).get('band')
            time_status=('通过' if time_control['magnitude_delta']<=2e-4 and time_control['phase_error_deg']<=1 and
                        (case=='reference' or time_control['relative_loss_change']<=.05) else '未通过') if time_control else '待完成'
            consistency=record.get('frequency_band_consistency') or {}
            magnitude=f'{100*consistency["magnitude_delta"]:.6g}' if consistency else '—'
            relative=f'{100*consistency["relative_loss_change"]:.4g}%' if consistency and case!='reference' else '不适用'
            lines.append(f'|{kind}／{case}|{time_status}|{magnitude}|{relative}|{record["window_delta_max"]:.6g}|'+('通过' if record['passed'] else '尚未通过')+'|')
    guided=root/'guided/base/L8/seed_101'
    if all((guided/c/'spectra.npz').exists() for c in ['reference','group']):
        group=dict(np.load(guided/'group/spectra.npz'));reference=dict(np.load(guided/'reference/spectra.npz'))
        if 'harmonic_band_T' in group and 'harmonic_band_T' in reference:
            actual=np.diagonal(group['T'],axis1=1,axis2=2)/np.diagonal(reference['T'],axis1=1,axis2=2)
            target=np.diagonal(group['harmonic_band_T'],axis1=1,axis2=2)/np.diagonal(reference['harmonic_band_T'],axis1=1,axis2=2)
            fig,axes=plt.subplots(2,2,figsize=(10,6),sharex=True)
            for i,mode in enumerate(['S0','A0']):
                f=group['frequencies']/1000
                axes[0,i].plot(f,(1-abs(target[:,i]))*100,'o-',label='独立频域计算')
                axes[0,i].plot(f,(1-abs(actual[:,i]))*100,'s--',label='当前瞬态计算')
                axes[0,i].set(title=mode,ylabel='同模态幅值下降（%）');axes[0,i].legend()
                relative=abs(abs(actual[:,i])-abs(target[:,i]))/np.maximum(np.maximum(abs(1-abs(actual[:,i])),abs(1-abs(target[:,i]))),1e-12)
                axes[1,i].plot(f,relative*100,'o-');axes[1,i].axhline(5,color='red',linestyle='--',label='5% 验收线')
                axes[1,i].set(xlabel='频率（kHz）',ylabel='损失相对差（%）');axes[1,i].legend()
                for ax in axes[:,i]:ax.grid(alpha=.2)
                lines += ['',f'{mode}：200 kHz 瞬态幅值下降 {100*(1-abs(actual[8,i])):.6g}%，频域 {100*(1-abs(target[8,i])):.6g}%；频带内最大损失相对差 {100*relative.max():.4g}%。']
            fig.suptitle('Lamb 波孔群：实际入射与无孔参考归一化；'+('导波全频带瞬态验收已通过' if guided_passed else '导波全频带瞬态验收进行中'))
            fig.tight_layout();fig.savefig(root/'latest_transient_comparison.png',dpi=160);plt.close(fig)
            lines += ['',f'![瞬态与频域对比]({(root/"latest_transient_comparison.png").resolve().as_posix()})']
    lines += ['', '边界选择已改为随模型尺寸调整容差，并核对实际覆盖长度；体波新增周期两侧相同网格的对照。旧漏选边界及未匹配周期网格结果保留在 diagnostics，不参与新结果拟合。',
        'A0／P 的时间步、完整观察窗和 17 个激励频带频点验证仍须以 validation.json 的验收记录为准，不能由此表的中心频率检查替代。',
        '求解加速试验未通过的版本已撤回；恢复为经连续求解对照验证的续算设置。',
        f'临时与恢复目录：E:/comsol_runtime/transient_matched/。C 盘当前可用 {payload["C_free_gib"]:.2f} GiB。','',
        f'![中心频率验证]({(root/"latest_center_validation.png").resolve().as_posix()})',
        f'![频带空间验证]({(root/"latest_band_validation.png").resolve().as_posix()})']
    sampling=root.parent/'transient_refined/guided_sampling/guided/sampling/L8/seed_101/reference'
    diagnostic=sampling/'sampling_checkpoint_diagnostic.json'
    reference_prechecks=root/'P_reference_prechecks.json'
    if reference_prechecks.exists():
        old=json.loads(reference_prechecks.read_text())
        lines += ['', '## P 波原 12 μs 参考窗的失败记录', '',
            f'原固定短窗的复信号窗口差为 {old["window_complex_delta"]:.6g}（门限 0.0002），尾部信号比例为 {old["tail_signal_fraction"]:.6g}（门限 0.0001），均未通过。',
            f'全频带瞬态／频域幅值差为 {100*old["band_harmonic_amplitude_delta"]:.6g} 个百分点，超过 0.02 个百分点目标。',
            '已保留原始记录，并改为从已有末态继续、独立检查参考窗稳定性；更新后的验收以正式 validation.json 为准。本段是历史失败证据，不代表延长后的结果。',
            f'[原短窗检查数据]({reference_prechecks.resolve().as_posix()})']
    full_sampling=root.parent/'transient_refined/guided_sampling/sampling_validation.json'
    if full_sampling.exists():
        sampling_result=json.loads(full_sampling.read_text())
        lines += ['', '## 完整接收采样对照', '',
            '同一求解轨迹、相同内部时间步与相同观察窗终点，对两种模态、匹配无孔／孔群进行完整分析。',
            '|输出间隔|相对频域的最大幅值差（百分点）|非零损失最大相对差|', '|---|---:|---:|']
        for row in sampling_result['rows']:
            lines += [f'|每 {row["output_stride"]} 个内部步|{row["magnitude_delta_pp"]:.6g}|{row["relative_loss_percent"]:.4g}%|']
        for check in sampling_result.get('one_vs_two_step_output_checks',[]):
            lines += ['',f'{check["mode"]} 每步／每两步输出的最大损失相对变化为 {100*check["relative_loss_change"]:.4g}%，'+('通过' if check['passed'] else '未通过')+'原有门限。']
        lines += ['', '保留每步输出。上述采样对照使用约 140 μs 轨迹；修正后的完整导波验收另使用 280 μs 并独立减半内部时间步。']
    elif diagnostic.exists():
        check=json.loads(diagnostic.read_text());common=sampling/'common_sample_check.json'
        lines += ['', '## 接收采样独立对照（尚在计算）', '',
            f'保持内部时间步不变，比较每步输出与每三步输出。目前已检查的无孔 S0 历史长度为 {check["record_end_us"]:.4g} μs。',
            f'三个检查频点的幅值差为 {", ".join(f"{v:.6g}" for v in check["magnitude_difference_percentage_points"])} 个百分点。',
            '这是同一已求解轨迹的输出采样诊断，不是完整双模态收敛验收；仍需完整时间窗及匹配孔群对照。']
        if common.exists():
            error=json.loads(common.read_text())['shared_receiver_history_relative_error']
            lines += [f'两次求解在共有采样时刻的接收历史相对差为 {error:.3g}，说明当前差异来自输出采样／积分，而非更改内部时间步。']
    (root/'LATEST_VALIDATION.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(payload,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
