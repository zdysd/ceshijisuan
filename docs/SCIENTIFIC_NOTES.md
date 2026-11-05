> 历史科学说明：其中旧机器路径与旧队列命令仅作追溯。新电脑请以根目录README、ENVIRONMENT和AGENTS为准。

# CFRP 随机场孔隙与超声传播

采用张翔（2011）的高斯自相关随机场、按孔隙率排序二值化方法。代码同时提供二维轮廓、三维随机孔隙、统计描述、厚度方向 SAFE 频散、COMSOL 6.2 局部模型、频域／瞬态求解和衰减验收。首轮计算记录见 `output/REPORT.md`；重新打开模型、同孔形单孔／孔群透射对照和修正后的验收见 `output/revalidation/REPORT.md`。

## 先区分示例与真实材料

用户尚无完整材料卡／失败的原模型。`config.json` 使用论文表2.1中的 Cxx=10.478 GPa、Cxy=3.029 GPa、Gxy=7.170 GPa、密度1815 kg/m³；其余刚度是明确给出的示例补全值，**不是实测数据**。COMSOL x,y 对应论文 x,z 截面，COMSOL z 为纤维延伸方向。完整矩阵按 xx,yy,zz,yz,xz,xy 的工程剪应变 Voigt 顺序输入，必须对称正定。

当前已执行的导波案例是均匀单向示例材料。SAFE 可以处理对称分层的三项平面应变刚度，存在面外剪切耦合／非对称铺层时会拒绝用 A0/S0 标签；一般铺层的 COMSOL 自动建模仍需材料卡、铺层和方向。二维截面、随机三维模型与真实CFRP不是自动等价的。

## 本机运行

科学计算依赖在项目 `.venv` 中，版本记录在 `requirements-lock.txt`。本机中文路径下 JPype 原生库不能可靠加载，COMSOL 脚本通过已有的 ASCII 路径 Python 环境启动，再添加项目依赖；没有修改那个环境。跨机器请设置 `CFRP_COMSOL_PYTHON`。

```powershell
./run.ps1 tests
./run.ps1 science generate --porosity 0.015 --distributions --output output/sample
./run.ps1 science generate --size 162 162 120 --correlation 20 16 12 --pixel 3 --output output/sample_3d
./run.ps1 science dispersion
./run.ps1 science calibrate --target-x 40 --target-y 20
./run.ps1 comsol --output output/bulk
./run.ps1 comsol --transient --output output/transient
./run.ps1 comsol --guided --frequency 200000 --output output/guided
./run.ps1 comsol --air --output output/air
./run.ps1 comsol --three-d --output output/three_d
./run.ps1 comsol --geometry output/sample/clustered/geometry.json --output output/clustered_bulk
./run.ps1 study validate --output output/validation
./run.ps1 study ensemble --samples 10 --output output/ensemble
./run.ps1 study screen --samples 10 --output output/screen
```

已保存的频域模型可不重算，仅重提取前行波及验收：

```powershell
& 'E:/MCP/venvs/comsol/Scripts/python.exe' extract_forward.py output/final_ensemble_5MHz
& 'E:/MCP/venvs/comsol/Scripts/python.exe' extract_guided.py output/comsol_guided_corrected
& 'E:/MCP/venvs/comsol/Scripts/python.exe' validate_forward.py output/final_validation
& 'E:/MCP/venvs/comsol/Scripts/python.exe' build_equivalent.py output/final_ensemble_5MHz
& 'E:/MCP/venvs/comsol/Scripts/python.exe' validate_holdout.py
& '.venv/Scripts/python.exe' make_report.py
```

`./run.ps1 study screen` 先做数值验收，通过后才运行0.1%、0.5%、1%、2%、3%、4%、5%的随机样本扫描；`make_report.py` 仅整理已保存的结果。完整批次耗时取决于孔的实际数量和轮廓。`run_study.py` 复用已完成的测量，但配置不同会拒绝复用。

## 数据与物理定义

- 所有孔隙率输入都是比例：`0.015` 代表1.5%，不是 `1.5`。尺寸、相关长度使用微米，传播坐标和衰减拟合使用米。
- 生成连续随机场时使用功率谱的平方根滤波实白噪声，并扩展区域后裁剪；固定随机种子。正比例改变标准差不会改变排名孔形。
- 二维矢量轮廓按实际面积校准，保留内部实体岛和边界截断孔。`pixel_pores.csv` 与 `vector_pores.csv` 不混用；它们对应不同几何离散。
- 三维通过三维随机场生成，细化同一连续场后提取封闭等值面，并按曲面体积校准。STL 曲面体积与体素计数不同；记录实际孔隙率及验收标志。有限ROI的闭合处理相当于保留薄的无孔外层，应做尺寸／分辨率检查。
- 标定的长度／宽度定义为未截断连通孔的水平／垂直包围盒尺寸。**未确认它们与论文的等效宽度 W_eq 完全同义**；演示标定验证算法，不声称复现其全部金相统计。
- 分布对照平移同一组不规则轮廓，不改面积、尺寸和方向，不允许重叠。四乘四分区的方差／均值比仅作聚集描述，不作普适显著性检验。
- 体波用COMSOL切线单元积分计算固相条件平均位移与截面平均法向应力，再分离前后行波。不要用固定数量的点跨微孔近似积分。导波联合位移与牵引应力分解A0/S0前后方向，保留非模态残差。模态衰减包括散射和转换，不能全部解释为吸收。
- COMSOL 使用 `exp(+iωt)`，SAFE 推导使用 `exp(+ikx-iωt)`，导波激励与投影已作共轭转换。

## 边界、网格与衰减

空孔以几何扣除构成，孔壁自由。空气对照单独使用压力声学与声固耦合，固体域不会使用零密度或零刚度。

三维默认使用MeshPy/TetGen生成包含全部孔面和孔腔的体网格，并作为COMSOL网格几何导入。它避免STL转CAD产生的细小面／薄缝；不删除孔隙。体网格总体积与目标孔腔体积必须一致。`--cad-three-d`仅保留作诊断，已知该路线可能产生退化网格。网格模型的几何实体属于`mesh1`，不应再按空的CAD序列检查域数。

体波参考阻抗取`Z=sqrt(rho*Cxx)`，COMSOL的`exp(+iωt)`约定下前行波满足`σxx=-iωZu`，故`u_forward=(u-σxx/(iωZ))/2`。这是参考介质投影，在多孔材料中仍需统计与尺寸验收。各向异性阻抗分量必须设置`GlobalSystem`；自动生成的`sys1`为边界坐标系，方向不同。

体波演示使用开放侧边界，只能作为有限区域诊断。`run_study.py` 采用上下匹配的横向周期单元、源端及接收端低反射条件，避免有限宽度产生板模态；周期性是**明确的边界建模假设**。两侧40 μm无孔缓冲带保证边界几何匹配，并补偿内部孔隙率以保持截面总孔隙率。必须完成宽度收敛后才能外推。

导波保留上下自由表面、板厚，端部用平滑吸收区。示例体波与导波使用了不同孔径统计，**不能直接把两套示例衰减数值作为相同孔隙的体波／导波对照**。真实对照应保持孔隙统计、传播方向和材料相同，再选取各自的频率／波长范围。

瞬态使用五周期 Hann 窗、Generalized-alpha、内部手动时间步、rhoinf=1、直接线性求解器。5 MHz 演示时间步是1/(60f0)，仅是窄带起点；应按实际激励的最高有效频率重新设置并减半检验。`gated_spectra` 保持绝对时间相位，接收时间窗必须涵盖完整首达脉冲并排除反射。

默认先保存全部接收历史NPZ，再仅保存最后一个全场及其导数，避免多GB模型文件。用`--keep-full-fields`保留全部全场。`compact_results.py`只允许压缩本项目`output/`内生成的模型；替换前重新打开并核对最后一帧。压缩存储不代表求解CPU时间自动减少。

频域衰减对至少三个截面拟合 `ln|A|=b-alpha*x`，先与无孔参考在相同位置比较。幅值衰减 Np/m 转 dB/m 乘8.685889638。保留负拟合值作为诊断，不截断成正数。随机集合先平均复信号，再求相干衰减；个体拟合均值与相干平均不是同一个估计量。

验收分别检查轮廓、网格、时间步、域尺寸、随机样本数。非零累计衰减变化的相对目标为5%；接近零时以0.01 Np为绝对门限。`exponential_fit_accepted` 另外检查拟合残差，`material_status` 独立标识材料可信度。**求解完成不等于得到可外推的衰减系数。**

## 等效化与适用范围

`homogeneous_propagation` 只是单一方向／模态的复波数指数传播函数，采用`exp(+ikx-iωt)`约定并拒绝负衰减；函数本身不执行域尺寸或随机集合验收。调用者必须先验收参数，COMSOL模型生成入口`build_equivalent.py`另行检查验收文件。它不会把体波标量衰减自动映射到 A0、S0。需要一般宏观有限元等效材料时，必须进一步辨识方向相关复刚度／非局部响应，并用未参与标定的模型验证。

`build_equivalent.py`核对频率和验收结果，再生成COMSOL一维复Helmholtz相干P波模型。默认验收文件来自`output/revalidation_corrected_domain/forward_gates.json`；旧验收遗漏实际孔隙率容差，已在复核中更正。示例仅验证1 mm传播，不能以近零绝对容差通过为由宣称厘米级衰减已收敛。其解析对照误差约0.064%，两组未见局部样本的平均幅值预测误差约0.40%；相位／复幅值误差更大，见报告。导波非模态残差仍较大，目前禁止将体波系数直接套用到A0/S0。

参考：张翔（2011）本地论文；COMSOL 6.2 安装文档及[瞬态波指南](https://www.comsol.com/support/knowledgebase/1118)、[Lamb波建模示例](https://doc.comsol.com/6.3/doc/com.comsol.help.models.aco.lamb_waves_ndt_plate/lamb_waves_ndt_plate.html)。

## 匹配孔隙的体波／Lamb 波瞬态入口

`run_transient.py` 保留原频域流程，新增二维 Solid Mechanics 瞬态对照。默认板厚／周期单元宽度 1 mm、孔隙段 8 mm、面积孔隙率 1.5%、约 80 μm 不规则孔径组、seed 101；体波 5 MHz，S0/A0 为 200 kHz。两类波读取相同孔形与位置，单孔取同一孔群中的最大孔。材料仍是示例补全参数，孔腔仍是二维面外延伸缺陷。

```powershell
./run.ps1 transient pilot
./run.ps1 transient spatial
./run.ps1 transient verify --no-frames
./run.ps1 transient report --render
# 仅在基础验收通过后执行；至少三个长度及三个完整种子族才能生成统计衰减。
./run.ps1 transient ensemble --samples 6 --no-frames
```

可用 `--kind bulk|guided`、`--case reference|single|group` 分批运行；`--output` 指定独立输出根目录。同目录缓存检查材料、网格、时间步、孔形及布局，改变配置需使用新的输出目录。默认输出为 `output/transient_matched/`，不覆盖早期算例。`spatial` 使用同一线性空间算子的中心频率谐波解，筛查网格减半及导波吸收区延长／体波周期宽度加倍；它不是整个脉冲频带的瞬态网格收敛证明。`verify` 先复用这些空间检查，再实际进行瞬态内部时间步减半。宽度变体保持孔径和孔隙率，并增加独立孔形样本，不能把变化完全归因于边界误差。

本流程启动 JVM 前把进程的 `TEMP`、`TMP`、JVM 临时目录、独立偏好目录及恢复目录指定到 `E:/comsol_runtime/transient_matched/` 下各自的进程子目录。运行时路径使用纯英文，避免本机 COMSOL 原生启动选项误解码中文路径；结果仍保存到项目 `output/`。启动时创建探测模型，核对其实际恢复目录，并在输出根目录的 `_runtime/session_<pid>.json` 保存路径记录。不写入用户全局 COMSOL 偏好设置。本次从 C 盘迁出的中断恢复文件在 E 盘保留，迁移执行了逐文件 SHA-256 校验。

五周期 Hann 窗对导波两个位移分量保留复模态的正交相位。最高有效频率取两种脉冲正交分量合计能量的 99.9% 频带上限；初始内部步长为 `1/(60*fmax)`，输出步长为 `1/(20*fmax)`。二者分别记录。脉冲中心模态不是全频带纯模态，接收谱在 0.7～1.3 倍中心频率内逐频求 SAFE 模态并分离前后行波；两次独立导波激励组成入射矩阵。四个接收变量为 `u,v,solid.sx,solid.sxy`，六个无孔截面完整历史保存为 `*_history.npz`。

同几何、同网格另做中心频率谐波求解，保存为 `*_harmonic.npz`，用于瞬态／频域交叉检查。频谱保留绝对时间；比较全记录与前 90% 记录，检查尾部位移能量和时间窗敏感性。导波报告模态拟合残差、远端反向入射和功率平衡；体波截面平均 P 通道不覆盖所有周期横向散射谐波，因此不以 P 通道能量缺口宣称材料吸收。

默认每次求解按最高有效频率每周期至少 8 帧保存真实瞬态位移网格 `*_fields.npz`，`report --render` 生成双分量 GIF 和接收波形 PNG。`*.mph` 保留最多 96 个代表性有限元快照，用于选定时刻的场检查；高频连续回放使用外部动画，避免稀疏快照的时间混叠。压缩前后核对所有保留时刻的完整截面场，`audit` 再重新打开模型复核。`--no-frames` 只保留最终有限元状态，但仍保留全部接收历史。会话已开启 E 盘自动恢复；`run.json` 记录当前模态与实际恢复目录，`*_solve_progress.json` 记录求解完成状态。

较早导波动画缓存仅含 96 帧，每个最高有效频率周期约 2.5 帧，满足 Nyquist 但回放较粗，仅供定性观察；新求解采用上述至少 8 帧的设置。`*_animation_sampling.json` 记录实际采样密度。可运行 `make_transient_previews.py` 将 GIF 转为逐帧核验的 MP4；本机视频工具位于 E:/comsol_runtime/media_tools，媒体临时文件也在 E 盘。

若尾部位移信号平方积分比例或时间窗敏感性检查失败，可用 `--duration-factor 2 --output output/transient_matched_longwindow` 延长记录；尾部比例不是总机械能。不能把短时间窗的谱当作已收敛稳态传递函数。基础验收失败时保留诊断结果，不生成材料衰减常数。

本机已用完整求解与分段续算对照验证重启状态衔接，见 `restart_check.json`。对已完成的体波单孔／孔群，可执行 `./run.ps1 transient extend --kind bulk --case group --duration-factor 2 --no-frames`，从保存的末端位移和速度继续计算尾波，不重新激励脉冲。原始结果保留，新记录位于 `bulk/window2/`；完整接收历史拼接后保留绝对时间，相应模型保留延长段末状态及其初始化解。延长窗结果复用已经衰减的匹配无孔参考，`window_validation.json` 单独报告窗口敏感性。

`summary.json` 与 `REPORT.md` 汇总体波 P、导波 S0/A0 的归一化透射，以及原始反射／转换通道。`validation.json` 明确区分尚未完成与已经失败的检查；基础数值验收未通过时拒绝运行 `ensemble`。`attenuation.json` 只在至少三个长度和三个完整种子族齐全时产生拟合与 Bootstrap 区间，单长度结果绝不伪造衰减系数。统计样本稳定性与网格收敛分开判断，有限样本拟合不自动标记为材料常数。

## A0 / P 收敛复验

`run_refined_pipeline.py all` 在 `output/transient_converged/` 续算。使用本机 COMSOL Python `E:/MCP/venvs/comsol/Scripts/python.exe` 启动；临时、自动恢复、队列日志都位于 E 盘。`LATEST_VALIDATION.md` 是带日期的检查快照，`pipeline_progress.json` 是当前阶段；目录名不代表验收已经通过。

新版按模型尺寸设置边界选择容差，核对完整边界长度，并对体波的周期两侧复制相同网格。空间检查覆盖 17 个激励频带频点；瞬态另检验内部时间步减半、完整波包、重启状态和频域一致性。续算关闭不必要的初始状态修正，但继续重新组装载荷和约束；未通过回归的常量载荷缓存试验不用于正式计算。

两类波的基础数值检查都通过后，程序才进入 8/16/24 mm 长度检查，再运行默认 12 个随机种子。长度扩展还会核对代表样本及拟合衰减斜率的网格敏感性。中心频率统计采用经瞬态交叉验证的频域算子；动画始终由真实瞬态位移帧生成，并直接导出逐帧核验的 MP4。统计结果单独写入 `attenuation_center_frequency.json` 和 `ATTENUATION_REPORT.md`，与单长度瞬态汇总区分。

任一数值验收失败，队列停止在失败阶段，不开始随机扩展。样本数稳定性、拟合适用性和材料参数可信度仍独立报告。体波宽度复制检查验证同一周期单元的数值一致性，不证明无限随机材料的代表体积已收敛。
