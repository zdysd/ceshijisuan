# Working agreement

- 当前仓库是可迁移版本；先阅读 README.md、docs/PROJECT_STATUS.md、docs/ENVIRONMENT.md。
- 默认只做离线检查。计算、补算、恢复旧队列均须有用户在目标电脑的明确指令。
- 不自动运行 resume_transient.ps1、restart_at_bulk_checkpoint.ps1、旧 P 波队列或任何 all/ensemble 扩展入口。
- 原冻结协议、几何、30个独立种子、阈值、窗口硬上限不得为通过验收而改变。
- 独立样本按样本、按 S0/A0 分别判断；单孔使用绝对幅值变化，孔群使用有符号损失。
- 复用前核对来源、端口、材料、签名、数据校验值；保留失败结果。
- 数值验收与实验 POD 分开。已有通过结果按兼容性复用。
- 源码变化会触发 execution_manifest 的引擎校验。不得简单覆盖旧哈希绕过检查；必须保留变更记录、兼容性依据并注册新版本。
- 不上传许可证、凭据、真实 .env、虚拟环境、生成结果和临时文件。未经授权不提交、推送、部署。
