# P1 剩余 7 个单孔 · 独立诊断

已得到逐项结论：6/7。状态：running。

48 mm，200 kHz。S0/A0 分别使用绝对同模态透射幅值变化判定，阈值为 0.1%。
在数值检查通过时，幅值变化减去数值敏感性仍大于 0.1% 判为超标，相加仍小于 0.1% 判为低于阈值。其余为待定。
本批为独立诊断。原 P1 失败结果保留，P2/P3 未启动。

| 算例 | 细化级别 | S0 幅值变化 ± 敏感性 | S0 判定 | A0 幅值变化 ± 敏感性 | A0 判定 |
|---|---:|---:|---|---:|---|
| `single_max_span_lower_r90` | 1 | 0.254410% ± 0.000695% | 超标 | 0.290811% ± 0.000396% | 超标 |
| `single_strongest_development_response_center_r0` | 1 | 0.004837% ± 0.000061% | 低于阈值 | 0.003802% ± 0.004177% | 低于阈值 |
| `single_strongest_development_response_upper_r0` | 1 | 0.025185% ± 0.000214% | 低于阈值 | 0.037019% ± 0.004739% | 低于阈值 |
| `single_strongest_development_response_lower_r0` | 1 | 0.024001% ± 0.000448% | 低于阈值 | 0.032876% ± 0.004497% | 低于阈值 |
| `single_strongest_development_response_center_r90` | 1 | 0.002653% ± 0.000025% | 低于阈值 | 0.003466% ± 0.000434% | 低于阈值 |
| `single_strongest_development_response_upper_r90` | 1 | 0.012210% ± 0.000079% | 低于阈值 | 0.015831% ± 0.002421% | 低于阈值 |
| `single_strongest_development_response_lower_r90` | — | — | 待计算 | — | 待计算 |

数据更新时间（UTC）：2026-11-09T10:57:43.562115+00:00。
