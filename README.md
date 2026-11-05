# CFRP 超声体波与 Lamb 波数值评估

本仓库用于 Windows 上的随机孔隙建模、体波及 S0/A0 导波评估。材料为示例补全参数，二维孔腔与理想端口；数值越阈能力不等于实验 POD。

## 在另一台电脑安装

1. 安装 Git、64 位 Python 3.12 和 COMSOL 6.2（计算时还需对应物理模块及有效许可证）。
2. 将私有仓库克隆到英文路径，例如 `D:/Projects/cfrp`，用本地 Codex 打开该目录。
3. 在项目目录的 PowerShell 执行 `./setup.ps1`。如果没有 Python Launcher，执行 `./setup.ps1 -Python C:/Python312/python.exe`。
4. 如需自定义临时盘，先设置 `$env:CFRP_RUNTIME_DIR = 'D:/comsol_runtime/cfrp'`。路径必须是英文，磁盘空间按待运行算例单独核对。
5. 执行 `./run.ps1 doctor` 和 `./run.ps1 tests`。两者均不启动 COMSOL 求解。

默认统一使用项目 `.venv/Scripts/python.exe`。如必须使用另一解释器，设置 `CFRP_COMSOL_PYTHON`，并在那个解释器中安装相同锁定依赖。不要复制旧电脑的 `.venv`。

## 迁移数据与当前进度

`data/portable_inputs/` 保存经过筛选的几何、冻结协议、既有端口数组和 P1 缓存。`tools/restore_inputs.py` 校验这些文件的 SHA-256，再恢复缺失的 `output/` 文件；不会覆盖目标电脑已有的结果。`output/` 本身仍不提交。

旧绝对路径由读取器映射到当前仓库，原协议、种子和几何文件的字节及哈希不改写。任何历史路径都优先解析为本仓库文件，不回读旧电脑目录。

来源为 2026-10-29 交接快照：P1 已有 14/22 个算例通过数值检查，剩余 8 个未完成；P2/P3 未开始。旧 P 波为 needs_adjustment，禁止自动重启。已通过数值检查不等于已达到检出阈值，seed102 仍低于主阈值。

本包未包含旧 P 波/瞬态的完整 MPH 场解与重启历史，不支持直接恢复旧 P 波队列。48 mm P1 通过签名与端口缓存复用；后续新算例必须在目标电脑完成原生求解预检后显式启动。

## 目标电脑上的计算入口

本次迁移准备不执行计算或补算。目标电脑先核对 COMSOL 6.2 build、许可、可用内存、核心数、磁盘与临时路径，再由用户明确授权一个最小原生求解。

授权后，48 mm 队列入口为 `./run.ps1 detectability run`。它继续未完成 P1，然后按 P1→P2→P3 门限执行；禁止修改 0.001 阈值、固定30种子或硬时间上限来凑通过。普通报告刷新可用 `./run.ps1 detectability report`。

注意：`detectability audit` 会重新生成部分审计/准备文件，并非纯只读检查；核验迁移输入请使用 doctor。

## 进一步说明

- [目标电脑恢复步骤与验收](docs/ENVIRONMENT.md)
- [实际阶段及待办](docs/PROJECT_STATUS.md)
- [上传范围及外部数据](docs/MIGRATION_MANIFEST.md)
- [原项目科学说明](docs/SCIENTIFIC_NOTES.md)
- [源码改动与缓存兼容依据](docs/CHANGELOG.md)

Word 报告生成依赖另装 `requirements-report.txt`。历史 PPT/视频/Word 渲染工具不在本次跨电脑验收范围。
