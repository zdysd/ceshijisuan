"""Maintain one result entry page while the authorized calculation runs."""
import ctypes,json,os,time,sys
from ctypes import wintypes
from datetime import datetime,timezone,timedelta
from pathlib import Path

ROOT=Path(__file__).resolve().parent/'output/transient_converged'
from cfrp_pores.portability import runtime_root
SCRATCH=runtime_root()/'queue_temp'
SCRATCH.mkdir(parents=True,exist_ok=True)
os.environ['TEMP']=str(SCRATCH);os.environ['TMP']=str(SCRATCH)


def read(name,default=None):
    path=ROOT/name
    if not path.exists():return default
    return json.loads(path.read_text(encoding='utf-8-sig'))


def alive(pid,expected_start=None):
    kernel=ctypes.windll.kernel32;kernel.OpenProcess.restype=ctypes.c_void_p
    kernel.CloseHandle.argtypes=[ctypes.c_void_p]
    kernel.GetExitCodeProcess.argtypes=[ctypes.c_void_p,ctypes.POINTER(wintypes.DWORD)]
    kernel.GetProcessTimes.argtypes=[ctypes.c_void_p]+[ctypes.POINTER(wintypes.FILETIME)]*4
    handle=kernel.OpenProcess(0x1000,0,int(pid))
    if not handle:return False
    try:
        code=wintypes.DWORD()
        if not kernel.GetExitCodeProcess(handle,ctypes.byref(code)) or code.value!=259:return False
        if expected_start:
            created,exited,cpu,user=[wintypes.FILETIME() for _ in range(4)]
            if not kernel.GetProcessTimes(handle,*[ctypes.byref(v) for v in [created,exited,cpu,user]]):return False
            actual=((created.dwHighDateTime<<32)|created.dwLowDateTime)/1e7-11644473600
            expected=datetime.fromisoformat(expected_start.replace('Z','+00:00')).timestamp()
            if abs(actual-expected)>10:return False
        return True
    finally:kernel.CloseHandle(handle)


def update():
    state=read('pipeline_progress.json',{});process=read('pipeline_process.json',{})
    stage=state.get('stage','starting');running=alive(process.get('launcher_pid',0),process.get('started'))
    sampling_process=read('sampling_process.json',{})
    sampling_running=alive(sampling_process.get('launcher_pid',0),sampling_process.get('started'))
    corrected=ROOT.parent/'transient_verified'
    corrected_process=read(Path('../transient_verified/continuation_process.json'),{})
    corrected_state=read(Path('../transient_verified/continuation_state.json'),{})
    corrected_running=alive(corrected_process.get('launcher_pid',0),corrected_process.get('started'))
    descriptions={'starting':'正在启动计算','waiting_for_existing_solves':'等待已有计算结束',
        'spatial_and_spectral_checks':'空间与频带验证','transient_validation':'瞬态收敛验证',
        'validation_complete':'基础数值验收结束','length_numerical_validation':'孔隙段长度与拟合斜率的网格验证',
        'random_seed_expansion':'随机样本扩展','complete':'队列已完成','needs_adjustment':'检查未通过，队列已停止扩展'}
    stamp=datetime.now(timezone(timedelta(hours=8))).isoformat(timespec='seconds')
    lines=['# 仿真结果入口','',f'更新时间：{stamp}','',f'**当前阶段：{descriptions.get(stage,stage)}。**']
    if not running and not corrected_running and stage not in ['complete','needs_adjustment']:
        lines+=['','计算进程已停止；已完成结果和 E 盘断点仍保留，需要续算。']
    if stage=='needs_adjustment':lines+=['',f'检查信息：{state.get("error","")}']
    if state.get('kind'):lines+=['',f'波型：{state["kind"]}；算例：{state.get("case", "—")}；变体：{state.get("variant", "—")}。']
    if state.get('length_mm') is not None:lines += [f'孔隙段长度：{state["length_mm"]:g} mm。']
    if state.get('seed') is not None:lines += [f'随机种子：{state["seed"]}；计划样本数：{state.get("samples", "—")}。']
    if sampling_running:
        for case,label in [('group','孔群'),('reference','无孔参考')]:
            folder=ROOT.parent/'transient_refined/guided_sampling/guided/sampling/L8/seed_101'/case
            path=folder/'run.json'
            if not path.exists():continue
            record=json.loads(path.read_text(encoding='utf-8'))
            if record.get('status')=='running':
                lines += ['',f'接收采样独立对照正在计算：{label}，{record.get("active_mode", "—")} 激励；内部步长保持不变，每步保存接收位移和应力。']
                break
    if corrected_state:
        corrected_labels={'waiting_for_sampling_controls':'等待完整接收采样对照',
            'checking_cached_guided_results':'核对已通过的导波参数，直接复用既有结果',
            'guided_native_output_and_long_window':'每步输出、280 μs 导波长窗及半时间步验证',
            'waiting_for_bulk_transient_controls':'导波已通过，等待 P 瞬态对照',
            'bulk_refinement':'P 波观察窗或时间步继续加密，既有结果已保留',
            'waiting_for_source_release':'P 算例已保存，等待原求解进程释放内存后继续验收',
            'numerical_acceptance_passed_starting_length_and_seed_controls':'全体验收通过，进行长度及随机样本扩展',
            'needs_adjustment':'仍有验收未通过，统计扩展被阻止','complete':'修正队列完成；样本稳定性仍以报告为准'}
        lines += ['',f'**修正后续算队列：{corrected_labels.get(corrected_state.get("stage"),corrected_state.get("stage"))}。**',
            f'[修正队列状态与新结果]({(corrected/"CONTINUATION_STATUS.md").as_posix()})',
            f'[额度刷新后的续算说明]({(ROOT.parent.parent/"CONTINUE_TRANSIENT.md").as_posix()})']
    if state.get('kind') in ['bulk','guided'] and state.get('variant') in ['base','time'] and state.get('case') in ['reference','group']:
        folder=Path(state['kind'])/state['variant']/'L8/seed_101'/state['case']
        record=read(folder/'run.json',{});mode=record.get('active_mode')
        if mode in ['P','S0','A0']:
            progress=read(folder/'blocks'/mode/'progress.json',{})
            lines += [f'当前激励：{mode}；已保存 {progress.get("completed",0)} 个连续时段。']
    gate=read('validation.json',{})
    lines += ['', '基础数值验收通过。' if gate.get('ready_for_ensemble') else '基础数值验收尚未全部通过；长度／随机样本扩展受验收门控制。',
        '临时文件、自动恢复和分段断点均使用 E:/comsol_runtime/transient_matched/ 或本结果目录。','']
    for name,label in [('LATEST_VALIDATION.md','最新验证结果'),('REPORT.md','已保存的波形、透射与反射报告'),
                       ('ATTENUATION_REPORT.md','长度／随机样本与衰减报告'),('attenuation_center_frequency.json','衰减数值及置信区间')]:
        path=ROOT/name
        if path.exists():lines.append(f'- [{label}]({path.as_posix()})')
    if not (ROOT/'ATTENUATION_REPORT.md').exists():lines+=['','衰减结果尚未生成。']
    videos=[]
    for kind,mode in [('bulk','P'),('guided','S0'),('guided','A0')]:
        path=ROOT/kind/'base/L8/seed_101/group'/f'{mode}_displacement.mp4'
        verification=read(path.relative_to(ROOT).with_name(f'{mode}_video_validation.json'),{})
        if path.exists() and verification.get('full_decode_verified'):
            label='数值验收进行中' if not gate.get('ready_for_ensemble') else '代表算例数值验收已通过'
            videos += ['',f'{mode} 孔群位移动画（{label}）：',f'![{mode} 位移动画]({path.as_posix()})']
    lines+=videos or ['','新版孔群位移动画尚未导出。']
    for name,label in [('latest_band_validation.png','频带验证'),('latest_transient_comparison.png','瞬态与频域对比'),('comparison_summary.png','对比结果'),('attenuation_length_fits.png','衰减长度拟合')]:
        path=ROOT/name
        if path.exists():lines+=['',f'![{label}]({path.as_posix()})']
    temporary=ROOT/f'RESULTS.{os.getpid()}.tmp'
    temporary.write_text('\n'.join(lines)+'\n',encoding='utf-8');temporary.replace(ROOT/'RESULTS.md')
    return (running and stage not in ['complete','needs_adjustment']) or sampling_running or corrected_running


if __name__=='__main__':
    sys.stdout.reconfigure(encoding='utf-8')
    while update():time.sleep(30)
