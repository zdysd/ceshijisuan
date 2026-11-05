"""Verify committed files and write a durable handoff after the workers stop."""
import argparse,hashlib,json,zipfile
from pathlib import Path
from datetime import datetime
import numpy as np

root=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--case',choices=['reference','group'],default='reference')
parser.add_argument('--variant',choices=['base','time'],default='base')
parser.add_argument('--minimum-completed',type=int,default=9)
args=parser.parse_args()
folder=root/f'output/transient_converged/bulk/{args.variant}/L8/seed_101/{args.case}'
blocks=folder/'blocks/P';p=json.loads((blocks/'progress.json').read_text(encoding='utf-8-sig'))
assert p['completed']>=args.minimum_completed,'Target block has not committed'
r=json.loads((folder/'run.json').read_text(encoding='utf-8-sig'));dt=r['signature']['layout']['output_dt']
end=0.;files=[]
for i in range(p['completed']):
    path=blocks/f'history_{i:03d}.npz'
    with np.load(path) as h:
        t=h['times'];assert np.isfinite(h['values']).all()
        np.testing.assert_allclose(t[0],end,rtol=1e-12,atol=1e-17)
        np.testing.assert_allclose(np.diff(t),dt,rtol=1e-9,atol=1e-17)
        end=float(t[-1])
    files.append(dict(path=str(path),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
assert end>0
checkpoint=blocks/p['checkpoint_file']
assert checkpoint.resolve().is_relative_to(blocks.resolve())
with zipfile.ZipFile(checkpoint) as z:assert z.testzip() is None
files.append(dict(path=str(checkpoint),bytes=checkpoint.stat().st_size,sha256=hashlib.sha256(checkpoint.read_bytes()).hexdigest()))
if args.case=='group':
    for i in range(p['completed']):
        path=blocks/f'P_{i:03d}_fields.npz'
        assert path.is_file(),f'Movie fields missing: {path}'
        with zipfile.ZipFile(path) as z:assert z.testzip() is None
        files.append(dict(path=str(path),bytes=path.stat().st_size,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
out=root/'output/transient_verified'
report=dict(stage='ready_for_shutdown',checked_at=datetime.now().astimezone().isoformat(),
    case=args.case,variant=args.variant,completed=p['completed'],end_time_us=end*1e6,checkpoint_file=str(checkpoint),files=files,
    scope='Complete saved receiver history and MPH archive integrity; numerical acceptance remains pending.',
    resume_command='powershell -NoProfile -File ./resume_transient.ps1 -ResumeAfterShutdown')
(out/'PAUSED_FOR_SHUTDOWN.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
text=f'''# 已安全停止，可正常关机

核验时间：{report['checked_at']}

P {args.variant} {args.case} 已完整保存至 {end*1e6:.9f} μs，断点 {checkpoint.name}。计算、等待和监视进程已停止；接收记录连续且有限，MPH 文件完整性核验通过。历史及断点都在 E 盘。关机准备完成不代表 P 数值验收通过。

开机后先检查真实进程，再在工作目录运行：

```powershell
powershell -NoProfile -File ./resume_transient.ps1 -ResumeAfterShutdown
```

参考窗口上限 50 μs、孔群 100 μs；上限内允许自动翻倍，触顶仍失败则输出诊断。参考基准已完成至 49.994 μs 并通过五项基准检查。本轮半时间步参考的全部时间块已保存，但后处理与频域复核须在恢复后完成；孔群半时间步尚未开始。P 全体验收仍待完成。A0／S0 已通过结果保留。全部通过后才做长度和随机样本扩展。48 mm 导波 P1 等待更细网格，须在恢复旧 P 队列后按资源策略接续。

下一次优先读取 PAUSED_FOR_SHUTDOWN.json、各进程元数据和 blocks/P/progress.json。run.json 中旧 running 状态属于停止前记录，不能据此判断进程存在。不要恢复旧自动重启守卫。
'''
(root/'SHUTDOWN_HANDOFF.md').write_text(text,encoding='utf-8')
(out/'shutdown_notification.txt').write_text(
    f'仿真已在 P {args.variant} {args.case} 的 {end*1e6:.6f} μs 完整断点停止。\n\n接收历史和模型文件已通过完整性检查，保存在 E 盘。\n现在可以正常关机了。\n\n下次开机后告诉我“已开机，继续”，将先检查进程和断点再续算。\n交接文件：{root / "SHUTDOWN_HANDOFF.md"}',encoding='utf-8')
with (root/'CONTINUE_TRANSIENT.md').open('r+',encoding='utf-8') as f:
    old=f.read();f.seek(0);f.write('> 最新状态：已按用户要求在完整断点停止，关机与开机续算优先按 [SHUTDOWN_HANDOFF.md](SHUTDOWN_HANDOFF.md) 执行。\n\n'+old);f.truncate()
print(json.dumps({k:v for k,v in report.items() if k!='files'}),flush=True)
