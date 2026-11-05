# 48 mm 导波后台计算交接（2026-10-24 23:20 CST 快照）

用户已授权启动此前等待资源的 48 mm 导波补算。旧 P 波任务保持原进程和断点，未重启或修改。新导波队列单会话串行、最多 2 核；仅调度策略改变，冻结的 200 kHz、S0/A0、48 mm、0.1% 示例阈值与 30 个独立种子不变。

- 新队列入口：`E:/MCP/venvs/comsol/Scripts/python.exe -u detectability_200kHz_L48.py run`。启动记录在 `worker_process.json`，本快照 launcher PID 31604、实际 Python PID 5700，实际 Python 已设 `BelowNormal` 优先级。队列正在等待空闲，尚未创建新的 COMSOL 会话。
- 已保存 P1 匹配无孔参考 `reference`、level 1、200 kHz S0 的原始端口数据 `solves/center/reference/level_1/base/200000_S0.npz`，401／801／1601 点均有限值；原生求解耗时 148.94 秒，见 `first_native_timing.json`。本算例和 P1 阶段都还没有完成双模态数值验收。续算应直接跳过已保存 S0，不重复求解。
- 首次原生会话的实际临时及模型恢复路径见 `_runtime/session_26456.json`，均已核验在 `E:/comsol_runtime/transient_matched/`；该会话现已关闭。其 COMSOL 进度记录在 `E:/comsol_runtime/transient_matched/solver_26456.log`。
- 旧 P 波 launcher PID 32492、实际 Python PID 9244，续算控制器 launcher PID 24816、实际 Python PID 8536；均为本快照，不应按历史 PID 直接操作。旧 P 波进度以 `output/transient_converged/pipeline_stdout.log` 和各块 `progress.json` 为准。
- 每次新求解前的调度规则见 `resource_policy.json`：连续 20 分钟没有用户输入、启动前至少 12 GiB 空闲物理内存、后续每次求解前至少 8 GiB。旧 P 波活跃时仅 P1 的网格级 0／1 可低优先级运行，而且每次求解前须有至少 12 GiB 空闲；更细的级 2／3 以及 P2/P3 等待旧队列完成。若用户重新操作电脑，当前求解完成并落盘后等待下一段空闲；不会中断正在求解的单个算例。
- 切换空闲调度时曾发生 Win32 优先级 API 参数类型错误，尚未启动新求解就退出；完整启动错误已保存在 `incidents/priority_startup_20261024_2315/`，修复后已实测 `BelowNormal`，不属于 P1 科学验收失败。
- 续接时先核对进程的创建时间、父子关系和命令行，再读 `worker_state.json`、`worker_error.json`、`first_native_timing.json`、`p1_validation.json`、`p2_validation.json`、`p3_validation.json`，检查已有 `solves/` 下完整签名与 `complete.json`，避免重复计算。
- 阶段依赖保持 P1→P2→P3；任一阶段 `needs_adjustment` 时保留失败证据，停止后续扩展。`REPORT.md` 与 `PLAN_AND_STATUS.md` 可用 `.venv/Scripts/python.exe detectability_200kHz_L48.py report` 刷新。
