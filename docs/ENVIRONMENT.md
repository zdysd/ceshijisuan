# 环境与目标电脑验收

## 支持范围

已准备并离线验证 Windows x64、Python 3.12。Python 科学依赖见 requirements-lock.txt；报告依赖为可选。COMSOL 6.2 的安装、build、Solid Mechanics 相关模块和许可证要在目标电脑独立核对。本仓库不分发 COMSOL 或任何许可证。

## 安装

克隆到英文路径，执行 setup.ps1：创建全新 .venv、安装锁定依赖、校验并恢复随仓库的紧凑输入、执行离线预检与单元测试。脚本任何一步失败即停止。

`CFRP_COMSOL_PYTHON`：可选的替代解释器；默认项目 .venv。替代环境须安装同一锁定依赖。

`CFRP_RUNTIME_DIR`：可选的英文绝对临时目录。默认 LOCALAPPDATA/cfrp_comsol_runtime；用户目录含非英文时必须设置到可写英文路径。求解时会校验实际恢复目录仍在该根目录内。

## 离线验收

- `./run.ps1 doctor`：检查 Python、平台、依赖、路径写权限、磁盘空闲量、冻结输入和引擎文件哈希。
- `./run.ps1 tests`：执行科学与工程回归测试，不启动 COMSOL。
- `./run.ps1 detectability report`：仅生成现有数据的报告。

以上通过仍不等于已验证许可证或实际 COMSOL 求解。doctor 的 JSON 明确返回 native_solver_tested=false 和 license_tested=false。

## 原生验收：仅在目标电脑授权后

先确认无重复求解进程；核对 COMSOL 6.2 build/物理模块/许可证及资源限制，使用独立输出目录做一个匹配无孔最小算例，确认模型启动、求解、保存、恢复路径和端口提取。记录命令、版本和结果，通过后才允许恢复48 mm队列。

不可直接恢复旧 P 波：本包没有旧 MPH 断点及完整历史。不得把新机器 PID 与原机历史 PID 混用。

48 mm P1 复用的是完整签名及校验过的端口缓存；P2/P3 未算。冻结协议保持原字节，旧根路径仅经 project_path 解析，不重写协议。

历史脚本中可能仍保留旧机器路径或文档渲染插件路径，标为存档工具；本次支持的迁移入口为 setup.ps1、run.ps1、doctor、tests 和 detectability 48 mm 主流程。
