# 执行范围与实际状态

更新时间：2026-11-04T07:52:41.325782+00:00

状态：imported_snapshot_not_running

旧P波任务的进程、窗口上限、断点和依赖关系保持原样。新COMSOL作业等待旧重计算释放资源后串行执行，最多2核；轻量审计和几何生成不启动JVM。

|阶段|范围|实际状态|
|---|---|---|
|P0|原始双入射矩阵重建、兼容性、12样本表、误差证据|已完成|
|P1|23个几何（含参考）、初始base/fine/absorber变体，401/801/1601端口提取|paused_snapshot|
|P2|至多4几何；17频点、网格/吸收区、双模态瞬态及半步；280→最多560 μs|pending|
|P3|冻结30种子，48 mm、200 kHz，至少成对网格，固定分母30|pending|
|P4|真实噪声、传感器、材料及实验清单|已列清单，实验待数据|

P1初始几何变体上限 69，双模态求解上限 138（未扣除可复用结果）。P3至少60个孔群网格变体、120个中心频率模态求解，匹配参考优先复用。开始各批次时另记录实测首例耗时，不将此前P波耗时当作导波FD耗时。

## P1模型与来源

|模型|类型|几何来源|
|---|---|---|
|reference|reference|匹配无孔矩形|
|group_seed_101|group|E:\Documents\ChatGPT\测试\output\guided_collective\d80\L48\seed_101\geometry\geometry.json|
|group_seed_102|group|E:\Documents\ChatGPT\测试\output\guided_collective\d80\L48\seed_102\geometry\geometry.json|
|group_seed_107|group|E:\Documents\ChatGPT\测试\output\guided_collective\extended_samples\d80\L48\seed_107\geometry\geometry.json|
|group_seed_109|group|E:\Documents\ChatGPT\测试\output\guided_collective\extended_samples\d80\L48\seed_109\geometry\geometry.json|
|single_max_area_center_r0|single|E:\Documents\ChatGPT\测试\output\detectability_200kHz_L48\single_controls\single_max_area_center_r0\geometry.json|
|single_max_area_upper_r0|single|E:\Documents\ChatGPT\测试\output\detectability_200kHz_L48\single_controls\single_max_area_upper_r0\geometry.json|
|single_max_area_lower_r0|single|E:\Documents\ChatGPT\测试\output\detectability_200kHz_L48\single_controls\single_max_area_lower_r0\geometry.json|
|single_max_area_center_r90|single|E:\Documents\ChatGPT\测试\output\detectability_200kHz_L48\single_controls\single_max_area_center_r90\geometry.json|
|single_max_area_upper_r90|single|E:\Documents\ChatGPT\测试\output\detectability_200kHz_L48\single_controls\single_max_area_upper_r90\geometry.json|
|single_max_area_lower_r90|single|E:\Documents\ChatGPT\测试\output\detectability_200kHz_L48\single_controls\single_max_area_lower_r90\geometry.json|
|single_max_span_center_r0|single|E:\Documents\ChatGPT\测试\output\detectability_200kHz_L48\single_controls\single_max_span_center_r0\geometry.json|
|single_max_span_upper_r0|single|E:\Documents\ChatGPT\测试\output\detectability_200kHz_L48\single_controls\single_max_span_upper_r0\geometry.json|
|single_max_span_lower_r0|single|E:\Documents\ChatGPT\测试\output\detectability_200kHz_L48\single_controls\single_max_span_lower_r0\geometry.json|
|single_max_span_center_r90|single|E:\Documents\ChatGPT\测试\output\detectability_200kHz_L48\single_controls\single_max_span_center_r90\geometry.json|
|single_max_span_upper_r90|single|E:\Documents\ChatGPT\测试\output\detectability_200kHz_L48\single_controls\single_max_span_upper_r90\geometry.json|
|single_max_span_lower_r90|single|E:\Documents\ChatGPT\测试\output\detectability_200kHz_L48\single_controls\single_max_span_lower_r90\geometry.json|
|single_strongest_development_response_center_r0|single|E:\Documents\ChatGPT\测试\output\detectability_200kHz_L48\single_controls\single_strongest_development_response_center_r0\geometry.json|
|single_strongest_development_response_upper_r0|single|E:\Documents\ChatGPT\测试\output\detectability_200kHz_L48\single_controls\single_strongest_development_response_upper_r0\geometry.json|
|single_strongest_development_response_lower_r0|single|E:\Documents\ChatGPT\测试\output\detectability_200kHz_L48\single_controls\single_strongest_development_response_lower_r0\geometry.json|
|single_strongest_development_response_center_r90|single|E:\Documents\ChatGPT\测试\output\detectability_200kHz_L48\single_controls\single_strongest_development_response_center_r90\geometry.json|
|single_strongest_development_response_upper_r90|single|E:\Documents\ChatGPT\测试\output\detectability_200kHz_L48\single_controls\single_strongest_development_response_upper_r90\geometry.json|
|single_strongest_development_response_lower_r90|single|E:\Documents\ChatGPT\测试\output\detectability_200kHz_L48\single_controls\single_strongest_development_response_lower_r90\geometry.json|

代表孔群选择：101, 102, 107, 109。中间值按排序后的较低中位秩取值，规则已冻结。
每个失败设置最多追加两级针对性细化；失败后写 needs_adjustment、保留全部版本并阻止依赖阶段。u_num为有限控制差异的保守估计，不是严格数学界。
运行入口：项目 .venv 执行 detectability_200kHz_L48.py audit；ASCII COMSOL Python 执行 detectability_200kHz_L48.py run。原脚本不被修改。
