# 可迁移版本变更

原工作区只读保留；此版本在独立待上传目录完善。

- 新增project_path：将源清单中的原机根目录映射到当前仓库，不改变冻结数据字节或协议内容。
- 新增CFRP_RUNTIME_DIR，解除48 mm主流程对E盘的固定要求；运行时仍检查ASCII路径及实际恢复位置。
- 分类器新增single参数；单孔按score=abs(loss)分类，孔群按有符号loss分类。原P1科学门、阈值与种子不变。
- NPZ改为临时写入、有限值校验、fsync、原子替换；复用前检查必要数组；新瞬态完成记录附历史SHA256。
- 新增setup、restore、doctor、回归测试、Git忽略规则及接手说明。
- 原execution_manifest保存在migration/original_execution_manifest.json；新引擎清单明确列出当前实现哈希及原始版本来源。

缓存兼容性：guided_collective.py未改，P1缓存的数值签名与NPZ原字节保留；新增代码只影响路径解析、报告分类和存储校验。输入和历史缓存校验通过后可在目标机复用。原生COMSOL跨机器一致性尚需目标机最小算例验证，不宣称已完成。

未处理：P波失败诊断、剩余P1/P2/P3及任何补算。
