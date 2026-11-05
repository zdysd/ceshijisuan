# 当前计算交接（2026-10-29）

核验时间：2026-10-29T15:03:09.377279+08:00

本次已在安全边界停止，当前可以关机。停止前没有旧 P 波进程，也没有 48 mm 导波的有限元求解在运行；48 mm worker 处于 `waiting_resources`，其打开的 COMSOL 客户端已在等待边界关闭。临时、模型恢复和结果路径继续使用 E 盘。

## 已保存结果

- 旧 P 波孔群半时间步已保存 22/22 个时间块，最后断点为 `restart_021.mph`。最终 P 数值门未通过，已达到参考 50 μs／孔群 100 μs 上限；长度和随机样本扩展没有启动。
- 48 mm 导波 P1 已保存 14 个通过的算例记录（4 个代表孔群、10 个单孔控制）。当前剩余 P1 level 3 算例未完成；P2 和固定 30 个独立样本没有启动。
- P1 的冻结参数仍为 200 kHz、S0/A0、48 mm、示例阈值 0.001。已完成缓存按签名复用，不能重复计算。

## 下次开机

先检查进程和本文件，再在工作目录运行：

```powershell
E:/MCP/venvs/comsol/Scripts/python.exe -u detectability_200kHz_L48.py run
```

它会读取已有 P1 缓存，从待完成的 level 3 继续，并遵守空闲和内存门限。完成后可刷新报告：

```powershell
.venv/Scripts/python.exe detectability_200kHz_L48.py report
```

旧 P 波已经进入 `needs_adjustment` 终态；不要自动重新启动旧 P 队列。若要重新处理其失败门，先保留当前证据并单独制定新的调整方案。

完整机器可读交接：[PAUSED_FOR_SHUTDOWN.json](output/detectability_200kHz_L48/PAUSED_FOR_SHUTDOWN.json)。
