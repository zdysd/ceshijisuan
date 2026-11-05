# 上传范围与迁移清单

## 建议上传

上传本可迁移目录中的源码、tests、tools、README、AGENTS、setup.ps1、run.ps1、requirements*.txt、config.json、docs、migration 与 data/portable_inputs。最终逐文件清单见 migration/upload_manifest.json。

data/portable_inputs 包含原 output/sample、guided_collective、detectability_200kHz_L48 的精选 JSON/NPZ/CSV/Markdown，以及旧体波终态验收摘要。它们约35MiB，是几何、协议、端口数据和P1缓存，不含COMSOL可执行程序、许可证或真实实验数据。

## 排除

.venv、output、__pycache__、.chart-data-*、.ppt-*、blobs、日志、锁文件、临时文件、MPH/MPHBIN、历史P波重启模型及完整瞬态全场均不上传。本次未删除原工作区中任何模型或断点。

## 恢复与完整性

运行 tools/restore_inputs.py 校验 migration/portable_data_manifest.json 中每个种子文件SHA256，再还原缺失文件到output。已有output文件不会被覆盖。doctor随后核验当前引擎、冻结协议、几何清单和独立样本几何。

.gitattributes 对冻结数据和Python源码使用 -text，防止Git自动换行转换破坏文件哈希。不要手工重新保存冻结JSON。

旧engine manifest保存在 migration/original_execution_manifest.json。可迁移版本另行登记源码哈希，保留原始版本关联；原协议哈希、几何字节与科学门限不变。

## 新电脑独立配置

COMSOL 6.2/build/物理模块/许可、内存和磁盘资源、CFRP_RUNTIME_DIR、可选的CFRP_COMSOL_PYTHON均需目标机核验。离线测试通过不等于原生求解通过。

## 发布状态

本目录是待上传版本，本次未初始化其Git历史、未暂存、未提交、未创建远程仓库或推送。
