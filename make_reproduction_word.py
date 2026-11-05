"""Chinese SOP and numerical reference data for the solution-free MPH package."""
from pathlib import Path
import json,csv,re
from collections import Counter
from docx import Document
from docx.shared import Inches,Pt,RGBColor,Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT,WD_CELL_VERTICAL_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

ROOT=Path(__file__).resolve().parent
DATA=ROOT/'output/guided_collective'
DEST=ROOT/'output/COMSOL_reproduction_20261014'
DOCX=DEST/'薄板孔隙群数值实验流程与复核手册.docx'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8'))


def font(run,size=11.5,bold=False,east='宋体'):
    run.font.name='Times New Roman';run.font.size=Pt(size);run.font.bold=bold;run.font.color.rgb=RGBColor(0,0,0)
    run._element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'),east)


def paragraph(text='',style=None):
    p=doc.add_paragraph(style=style);font(p.add_run(text))
    if text.startswith('表'):p.paragraph_format.keep_with_next=True
    return p


def heading(text,level=1):
    p=doc.add_heading(text,level)
    for r in p.runs:font(r,16 if level==1 else 13,True,'黑体')
    return p


def table(headers,rows,widths):
    t=doc.add_table(rows=1,cols=len(headers));t.alignment=WD_TABLE_ALIGNMENT.CENTER;t.autofit=False
    for col,w in zip(t.columns,widths):col.width=Cm(w)
    def fill(cell,text,header=False,index=0,col=0):
        cell.width=Cm(widths[col]);cell.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
        p=cell.paragraphs[0];p.paragraph_format.space_after=Pt(2);p.paragraph_format.space_before=Pt(2)
        if header:p.paragraph_format.keep_with_next=True
        p.paragraph_format.line_spacing=1.08;p.alignment=WD_ALIGN_PARAGRAPH.LEFT if col==0 else WD_ALIGN_PARAGRAPH.CENTER
        run=p.add_run(str(text));font(run,10,header)
        if header:run.font.color.rgb=RGBColor(255,255,255)
        pr=cell._tc.get_or_add_tcPr();sh=OxmlElement('w:shd');sh.set(qn('w:fill'),'234E70' if header else ('F1F5F8' if index%2 else 'FFFFFF'));pr.append(sh)
        margins=OxmlElement('w:tcMar')
        for side in ['top','left','bottom','right']:
            el=OxmlElement('w:'+side);el.set(qn('w:w'),'80');el.set(qn('w:type'),'dxa');margins.append(el)
        pr.append(margins)
    for i,h in enumerate(headers):fill(t.rows[0].cells[i],h,True,col=i)
    repeat=OxmlElement('w:tblHeader');t.rows[0]._tr.get_or_add_trPr().append(repeat)
    for index,row in enumerate(rows):
        cells=t.add_row().cells
        for col,value in enumerate(row):fill(cells[col],value,index=index,col=col)
        keep=OxmlElement('w:cantSplit');t.rows[-1]._tr.get_or_add_trPr().append(keep)
    borders=OxmlElement('w:tblBorders')
    for side in ['top','left','bottom','right','insideH','insideV']:
        b=OxmlElement('w:'+side);b.set(qn('w:val'),'single');b.set(qn('w:sz'),'4');b.set(qn('w:color'),'D9D9D9');borders.append(b)
    t._tbl.tblPr.append(borders)
    paragraph().paragraph_format.space_after=Pt(2)
    return t


def mr(text):
    r=OxmlElement('m:r');t=OxmlElement('m:t');t.text=text;r.append(t);return r


def fraction(top,bottom):
    f=OxmlElement('m:f')
    for tag,value in [('num',top),('den',bottom)]:
        n=OxmlElement('m:'+tag);n.append(mr(value));f.append(n)
    return f


def power(base,exponent):
    if not isinstance(base,str):
        delimiter=OxmlElement('m:d');pr=OxmlElement('m:dPr')
        for tag,val in [('begChr','('),('endChr',')')]:
            a=OxmlElement('m:'+tag);a.set(qn('m:val'),val);pr.append(a)
        delimiter.append(pr);inside=OxmlElement('m:e');inside.append(base);delimiter.append(inside);base=delimiter
    s=OxmlElement('m:sSup');e=OxmlElement('m:e');e.append(base if not isinstance(base,str) else mr(base));s.append(e)
    p=OxmlElement('m:sup');p.append(mr(str(exponent)));s.append(p);return s


def equation(*parts):
    p=doc.add_paragraph();p.alignment=WD_ALIGN_PARAGRAPH.CENTER
    m=OxmlElement('m:oMath')
    for part in parts:m.append(mr(part) if isinstance(part,str) else part)
    p._p.append(m);return p


def figure(name,caption):
    p=doc.add_paragraph();p.paragraph_format.keep_with_next=True
    p.add_run().add_picture(str(DATA/name),width=Cm(17.0))
    c=paragraph(caption);c.alignment=WD_ALIGN_PARAGRAPH.CENTER
    for r in c.runs:font(r,10)


def fmt(v,n=5):return f'{v:.{n}f}'
def pct(v,n=3):return f'{v*100:.{n}f}%'


DEST.mkdir(parents=True,exist_ok=True)
inventory=read(DEST/'source_inventory.json');main=[r for r in inventory if r['category']=='main_models']
cfg=read(ROOT/'config.json');settings=read(DATA/'study_settings.json')
stats=read(DATA/'extended_attenuation.json');checks=read(DATA/'numerical_checks.json');sweep=read(DATA/'sweep_summary.json')
single=read(DATA/'single_controls_summary.json');threshold=read(DATA/'threshold_sensitivity.json')
doc=Document();section=doc.sections[0]
for border in list(doc.styles.element.xpath('.//w:pBdr')):border.getparent().remove(border)
section.page_width=Inches(8.5);section.page_height=Inches(11)
section.top_margin=Cm(1.9);section.bottom_margin=Cm(1.9);section.left_margin=Cm(2.0);section.right_margin=Cm(2.0)
section.header_distance=Cm(.8);section.footer_distance=Cm(.8)
normal=doc.styles['Normal'];normal.font.name='Times New Roman';normal.font.size=Pt(11.5)
normal._element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'),'宋体')
normal.paragraph_format.line_spacing=1.18;normal.paragraph_format.space_after=Pt(6)
for name in ['Title','Subtitle','Heading 1','Heading 2','Heading 3']:
    st=doc.styles[name];st.font.color.rgb=RGBColor(0,0,0)
    st._element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'),'黑体')
head=section.header.paragraphs[0];font(head.add_run('薄板孔隙群数值实验  COMSOL 复核手册'),9)
foot=section.footer.paragraphs[0];foot.alignment=WD_ALIGN_PARAGRAPH.CENTER;font(foot.add_run('第 '),9)
field=OxmlElement('w:fldSimple');field.set(qn('w:instr'),'PAGE');foot._p.append(field);font(foot.add_run(' 页'),9)
title=doc.add_paragraph('薄板孔隙群数值实验流程与复核手册','Title')
for r in title.runs:font(r,23,True,'黑体')
paragraph('COMSOL 6.2  二维 Lamb 波  无求解结果模型包','Subtitle')
paragraph('编制日期 2026年10月14日')
paragraph('本手册用于在你自己的电脑上重新计算孔隙群模型，并与已有数值结果逐项核对。这里的“实验”是数值实验，材料为明确标注的示例 CFRP 参数，不是实验室实测材料。交付的 MPH 文件保留模型、参数、网格和求解器设置，不保留场解；参考数值保存在本手册和小型数据表中。')
paragraph('当前结果支持在指定透射幅值判据下观察到孔群累积响应，但尚不支持单孔在所有通道均不可检测，也未达到衰减系数随随机样本数变化小于5%的目标。复核时应同时检查通过项和未通过项。')
heading('1 复核目标与操作顺序')
for text in [
    '确认模型版本与材料假设。使用 COMSOL 6.2 或可兼容读取该版本的环境；本轮主模型采用固体力学、二维平面应变和频域求解。',
    '先计算无孔参考，再计算最小、中位和最大单孔；采用同一频率、入射模态、网格级别和接收位置。',
    '计算20、40、80微米孔径组，分别进行50、100、200千赫激励；各基本孔群使用6个随机样本。',
    '对80微米组计算8、24、48毫米孔隙段；200千赫使用12组随机种子，分别保存S0和A0激励的数据。',
    '在孔前与孔后无孔截面分离前后行波，求入射归一化的透射与反射矩阵，再扣除无孔参考。',
    '拟合孔隙段长度与相干幅值，计算衰减系数及95%统计区间，检查网格、样本数和阈值敏感性。']:
    paragraph(text,'List Number')

heading('2 材料坐标与几何参数')
paragraph('COMSOL 的 x 为传播方向，y 为板厚方向，z 为面外方向。按原方法验证时，面外方向对应纤维方向；这不是一般铺层板的完整三维验证。二维孔腔在物理意义上沿面外延伸。应变采用工程剪切记法，刚度排列为 xx、yy、zz、yz、xz、xy。')
table(['项目','数值或设定'],[
    ['密度','1815 kg/m³'],['材料模型','线弹性各向异性  平面应变'],['本底材料阻尼','0  仅吸收区附加阻尼'],
    ['板厚','1.000 mm'],['总长度',f'{settings["length"]*1000:.9f} mm'],['孔隙段起点',f'{settings["start"]*1000:.9f} mm'],
    ['吸收区起点',f'{settings["sponge"]*1000:.9f} mm'],['吸收区长度',f'{settings["absorb_length"]*1000:.9f} mm'],
    ['孔隙段长度','8、24、48 mm'],['孔隙段面积孔隙率','1.5%  不等于整块338 mm板的孔隙率'],['孔壁与上下表面','自由边界  不赋零密度固体'],
    ['孔隙内部','几何扣除形成空腔  主模型不求解空气'],['孔至板表面及孔隙段端部缓冲','生成孔群时50 μm']], [6,11.5])
paragraph('表1  示例刚度矩阵  单位 GPa  行列顺序均为 xx、yy、zz、yz、xz、xy')
labels=['xx','yy','zz','yz','xz','xy'];matrix=cfg['stiffness_GPa_voigt_xx_yy_zz_yz_xz_xy']
table(['分量']+labels,[[label]+[f'{v:g}' for v in row] for label,row in zip(labels,matrix)],[2.1]+[2.5]*6)
paragraph('参数来源说明：配置文件仅将 Cxx、Cxy、Gxy 和密度标为取自张翔论文表2.1；其他矩阵元素是示例补全假设。替换真实材料或铺层后，必须重新计算模态形状和波数，不能仅修改材料卡后继续沿用本包内置的激励剖面。刚度实际设置在固体力学的线弹性材料节点中。')
heading('3 孔隙生成与孔径分组')
paragraph('基准方法为高斯型自相关函数生成随机场，按数值排序选取最低的指定比例进行二值化，再提取不规则孔轮廓。相关长度 a、b 不是单个孔的长短半轴；固定孔隙率时正比例缩放随机场不改变排序。')
equation('R(Δx,Δz)=exp[−',power(fraction('Δx','a'),2),'−',power(fraction('Δz','b'),2),']')
table(['基准生成参数','设置'],[['相关长度 a、b','20 μm、10 μm'],['基准像素尺寸','2 μm'],['基准随机场区域','1000 μm × 1000 μm'],['方向','0°'],['轮廓库','取已生成不规则轮廓中面积最大的64个'],['直径定义','等效面积直径  不是长轴或费雷特直径']],[6,11.5])
equation('d = ',power(fraction('4A','π'),'1/2'))
paragraph('本轮受控分组保留轮廓形状，按名义等效直径的约0.8～1.2倍抽样并缩放，随后进行总面积归一化，使孔隙段孔隙率准确达到1.5%。孔轮廓重新无重叠放置，并保留50微米缓冲层。因此，本轮最终位置相关性不等于原高斯随机场的自相关函数。不能将本轮结果解释为纯粹只改变孔径而其他统计量完全不变。')
table(['名义组','实际直径范围','8 mm孔数','随机种子'],[['20 μm','15.97～24.06 μm','377','101～106'],['40 μm','32.01～48.63 μm','94','101～106'],['80 μm','63.03～99.52 μm','24','101～112用于扩展']], [3,5.2,3.2,6.1])
paragraph('80微米组的24和48毫米孔隙段分别含71和141个孔。不同长度使用同一种子重新生成，并非截取或延长同一块孔隙试件。单孔对照包括代表性最小、中位、最大孔，靠近上下表面的最大孔，以及从全部组样本中选取的实际最大等效直径孔。单孔三位置检查不能穷尽所有形貌与位置。')

heading('4 激励边界网格和求解器')
table(['设置','数值或操作'],[['激励形式','频域谐波位移  非时域脉冲'],['频率','50、100、200 kHz'],['入射模态','S0、A0分别激励'],['位移幅值','u_amp = 1×10⁻⁸ m  即10 nm'],['左端边界','按对应厚度模态剖面施加复位移'],['右端边界','低反射边界  前置渐变吸收区'],['标准最大网格尺寸','150 μm'],['加密最大网格尺寸','75 μm'],['标准最小网格尺寸','min 最大网格尺寸/30 与最小等效孔径/30'],['孔边加密检查','最小孔径/60  与背景75 μm配合'],['网格增长率与曲率系数','1.35 与0.35'],['求解器','线性频域  完全耦合  PARDISO直接求解器'],['时间步','本轮频域计算不设置时间步']], [6,11.5])
paragraph('这里的最小网格尺寸是网格控制参数，不表示所有孔边单元都严格等长。几何轮廓分辨率与有限元网格分开控制。已有网格随文件保留，复核原始数值时优先使用该网格。')
paragraph('吸收区采用三次渐变质量比例阻尼。孔隙区内无本底阻尼。下式中 x_s 为吸收区起点，L_a 为吸收区长度，f 为所选激励频率；x不超过x_s时该项为0，刚度比例阻尼为0。')
equation('αₘ(x)=6πf ',power(fraction('x−xₛ','Lₐ'),3),'  （x>xₛ）')
paragraph('默认相对几何修复容差曾合并邻近孔。已修正的20微米seed106和40微米seed105，以及部分扩展样本采用1 nm绝对修复容差重新构建。导出的本轮主模型使用修正后的版本；对应两个旧版本归入已知问题历史目录。')
table(['频率','S0波长','A0波长'],[['50 kHz','45.999 mm','8.918 mm'],['100 kHz','22.995 mm','6.163 mm'],['200 kHz','11.488 mm','4.175 mm']],[4.5,6.5,6.5])
paragraph('参考模态由120个厚度单元的半解析有限元模型得到；240单元加密后最大波长变化约0.000275%。第一高阶模态截止频率约0.994 MHz。孔径远小于波长不能直接等价于不可检测。')

heading('5 在 COMSOL 中重新计算')
numbering=doc.part.numbering_part.element
num_id=max(int(x.get(qn('w:numId'))) for x in numbering.findall(qn('w:num')))+1
num=OxmlElement('w:num');num.set(qn('w:numId'),str(num_id))
abstract=OxmlElement('w:abstractNumId');abstract.set(qn('w:val'),'7');num.append(abstract)
override=OxmlElement('w:lvlOverride');override.set(qn('w:ilvl'),'0');start=OxmlElement('w:startOverride');start.set(qn('w:val'),'1');override.append(start);num.append(override);numbering.append(num)
for text in [
    '将整个模型包解压至本地目录，建议路径简短。优先从 main_models/reference/model.mph 开始；该文件是本轮主参考，打开时没有求解结果。先另存一份工作副本。',
    '在“全局定义 → 参数”中设置 mode_id。0表示S0，1表示A0。设置 freq_id 选择该文件允许的频率。f0由公式自动生成，不要用固定数值覆盖该公式。',
    '检查该文件的频率编号。三频文件中1、2、3分别为50、100、200 kHz；双频文件中1、2为100、200 kHz；仅200 kHz的文件只有编号1。每个文件的参数说明及同名 parameters.json 均列出映射。其他编号不受支持。',
    '检查“固体力学”的指定位移节点（内部标记source），源表达式根据频率编号与模态编号调用内嵌插值函数。请保留这些函数。仅改变频率时源形状也随编号同步切换。',
    '检查网格与“研究 → 频域”的频率列表 f0，点击“计算”。计算完成后，分别导出孔前、孔后截面的复位移与应力。',
    '把 mode_id 改为另一个值，重新计算并导出第二次激励的数据。两次独立入射都需要保留；单次A0峰值不能替代完整模态矩阵。',
    '对孔群与对应无孔参考执行完全相同的操作。24、48毫米长度分别使用 reference_L24、reference_L48；加密算例使用 checks 目录中的相应加密参考。',
    '若进行参数化扫描，可将 mode_id 的0、1和本文件有效 freq_id 编号设为全组合扫描。首次复核建议逐项计算，先对照单项数据再批量运行。']:
    p=paragraph(text,'List Number')
    np=p._p.get_or_add_pPr().get_or_add_numPr()
    np.get_or_add_ilvl().val=0;np.get_or_add_numId().val=num_id
paragraph('几何中的孔位置和轮廓已固定。随机种子、孔径组、孔隙率和孔数记录在同名参数说明中；修改文字记录不会重新生成孔隙。要改变这些几何统计量，需要修改轮廓或重新生成几何。历史模型保留原来的激励和参数，不自动获得本轮频率与模态切换功能。')
heading('6 接收截面与模态分离')
table(['孔隙段','孔前x坐标 mm','孔后x坐标 mm'],[
    ['8 mm','79.998904\n81.998904\n83.998904','115.998904\n117.998904\n119.998904'],
    ['24 mm','同上','123.998904\n125.998904\n127.998904'],
    ['48 mm','同上','147.998904\n149.998904\n151.998904']], [3,7.25,7.25])
paragraph('每个截面沿y方向从1×10⁻⁸ m到0.001−1×10⁻⁸ m均匀取401点。导出 u、v、solid.sx、solid.sxy 的复数值，或分别导出实部和虚部共8列；保留坐标、频率和入射模态。COMSOL中可使用“结果 → 数据集 → 二维截线”，再导出数据；网格采样点数与这里的401点不同，应另做统一重采样。')
paragraph('将S0前行、S0后行、A0前行、A0后行的参考位移与牵引应力组成模态基。应力以2πf√(ρCxx)缩放，联合最小二乘求四个复系数。截面均位于无孔区；先按各模态波数补偿传播相位，再对三个截面平均。COMSOL采用 exp(iωt) 约定，本包源剖面与前后行基采用匹配的复共轭约定。')
equation('a₊ = ',fraction('1','N'),'∑ a₊(xⱼ) exp(ikxⱼ)')
paragraph('两次激励的实际孔前前行系数组成入射矩阵 Fin；孔后前行系数为 Fout，孔前后行系数为 Bpre。透射矩阵与反射矩阵如下。该归一化能够计入实际混合入射和源端反馈；不能先把不同位置的复数幅值不加相位校正直接平均。')
equation('T=Fout ',power('Fin','−1'),'     R=Bpre ',power('Fin','−1'))
paragraph('同一模态的T对角元除以匹配无孔参考的T对角元，得到归一化复透射 t。主观测量仅为幅值变化；反射、相位与转换模态分别报告，不以该幅值判据推断所有通道的可检测性。')
equation('δsingle=|1−|t||     Δgroup=1−|mean(t)|')
heading('7 衰减与统计计算')
paragraph('孔隙段长度取8、24、48毫米，拟合复透射集合平均的模值随长度的变化，并保留截距，以减少界面项混入斜率。接收截面位于孔外无孔段，不能把这些孔外截面的空间斜率当作孔隙区衰减。')
equation('ln|mean(t(L))|=b−αL')
equation('αdB/m=',fraction('20','ln 10'),' αNp/m     αdB/mm=',fraction('αdB/m','1000'))
paragraph('同时计算逐样本幅值定义的衰减，即先对 ln|t| 作集合平均再拟合，区分样本之间的相位相消与单个试件的幅值损失。对随机种子族进行1000次Bootstrap重采样，报告均值、标准差及95%区间。区间仅反映本轮样本统计变化，不包含材料误差和全部数值误差。')
paragraph('用户选择的主阈值为0.1%，另检查0.01%和1%。判断单孔时以已测单孔中的最大幅值变化为依据，并参考加密差；判断孔群时比较统计下界及数值差。代表性网格差不是每个几何和位置的严格误差上界，尤其不能把0.01%级别判断当作已充分验证。')

heading('8 已有数值数据供复核')
paragraph('表2  8毫米孔隙段的相干集合幅值下降  每组6个随机样本  括号内为95%统计区间')
rows=[]
for d in [20,40,80]:
    for f in [50000.,100000.,200000.]:
        pair=[next(r for r in sweep if r['diameter_group_um']==d and r['frequency_hz']==f and r['mode']==m) for m in ['S0','A0']]
        cells=[]
        for r in pair:cells.append(pct(r['amplitude_reduction'])+'\n('+ '～'.join(pct(v) for v in r['ci95_amplitude_reduction'])+')')
        rows.append([f'{d} μm',f'{f/1000:g} kHz']+cells)
table(['孔径组','频率','S0幅值下降','A0幅值下降'],rows,[2.5,2.5,6.25,6.25])
figure('孔径频率扫描.png','图1  孔径与频率扫描  负下降值原样保留')
paragraph('表3  80微米组200千赫的长度对照  12组随机样本')
table(['指标','S0','A0'],[
    ['48 mm幅值下降']+[pct(r['longest_segment_contrast']['amplitude_reduction']) for r in stats],
    ['上述下降的95%区间']+['～'.join(pct(v) for v in r['longest_segment_contrast']['ci95_amplitude_reduction']) for r in stats],
    ['α  Np/m']+[fmt(r['fit']['alpha_np_m']) for r in stats],
    ['α的95%区间  Np/m']+['～'.join(fmt(v) for v in r['ci95_alpha_np_m']) for r in stats],
    ['逐样本幅值定义α  Np/m']+[fmt(r['typical_fit']['alpha_np_m']) for r in stats],
    ['三长度拟合R²']+[fmt(r['fit']['r_squared'],4) for r in stats],
    ['样本数6→12时α变化']+[pct(r['sample_count_relative_change'],1) for r in stats],
    ['样本数变化≤5%','未通过','未通过']], [6.5,5.5,5.5])
figure('孔群长度与统计区间.png','图2  12组样本的长度对照  灰线连接相同种子的独立长度模型')
paragraph('表4  单孔对照中的最大幅值变化  百分比  含最小中位最大和移动位置检查')
rows=[]
for d in [20,40,80]:
    for f in [50000.,100000.,200000.]:
        pair=[next(r for r in single if r['diameter_group_um']==d and r['frequency_hz']==f and r['mode']==m) for m in ['S0','A0']]
        rows.append([f'{d} μm',f'{f/1000:g}']+[pct(r['max_transmission_change'],6) for r in pair]+[pct(r['max_reflected_modal_change'],6) for r in pair])
table(['孔径组','kHz','S0透射','A0透射','S0反射','A0反射'],rows,[2.5,2,3.25,3.25,3.25,3.25])
figure('单孔透射与反射判据.png','图3  透射与反射判据需分别判断')
paragraph('表5  80微米组200千赫48毫米孔隙段的阈值敏感性')
table(['模态','阈值','已测单孔透射低于阈值','孔群下界扣除代表网格差后超阈值','单孔反射低于阈值'],[
    [r['mode'],pct(r['threshold'],2)]+['是' if r[k] else '否' for k in ['single_below','group_ci_lower_minus_representative_mesh_above','reflected_single_below']] for r in threshold],[2,2,4.1,5.3,4.1])

heading('9 验证状态与复核判据')
table(['检查','已有结果','含义'],[['模态相速度交叉检查','最大相对差0.00161%','当前示例板频段内一致'],['互易与最小二乘投影','幅值差最大0.00658个百分点','代表性A0场的独立核对'],['无孔端口场拟合残差','全记录最大0.459%','不是接收幅值误差'],['模态功率平衡','最大偏差0.1024%','仍需独立数值误差判断'],['随机样本数稳定性','S0约43%  A0约10%','均未达到5%目标'],['单孔全通道不可检测','不成立','部分反射超过0.1%']], [5.5,6,6])
table(['频率','模态','代表性α网格变化','5%目标'],[[f'{r["frequency_hz"]/1000:g} kHz',r['mode'],pct(r['relative_alpha_change'],2),'通过' if r['relative_alpha_mesh_passed'] else '未通过'] for r in checks['alpha_checks']],[4,3,6,4.5])
paragraph('上述网格对照仅使用seed101，不能替代整个随机集合的网格收敛。轮廓逐级加密、独立吸收区长度变化、三维薄板、真实铺层以及实测噪声下的检测性能尚未完成验证。当前α只能作为本次参数和长度范围内的局部估计，不能作为已收敛材料常数。')

heading('10 文件结构与核对方法')
table(['目录或文件','用途'],[['main_models','本轮84个薄板主模型  每个文件配同名parameters.json'],['historical_models','早期体波瞬态三维及其他验证版本  保留原设置'],['historical_invalid','两例已知孔隙合并的旧版本  仅用于追溯'],['模型清单.csv','全部模型与原始文件映射  导出状态和文件大小'],['参考数据','本轮汇总数值  不包含位移场求解向量'],['export_main.json 与 export_history.json','逐模型参数网格数量及重新打开后的空解检查记录'],['recompute_smoke_check.json','对导出副本重新求解的抽样一致性核对  不回写场解']], [6,11.5])
paragraph('主模型的几何与激励数据已内嵌，复算不依赖原电脑的孔轮廓路径。部分历史三维导入模型带有assets目录；请连同所在文件夹一起复制。若重建导入几何或网格时软件提示找不到输入文件，在相应导入节点重新选择同目录assets中的文件。历史模型未自动修正早期物理假设或已知局限，不能与本轮通过核对的主模型混用。')
paragraph('使用所有主模型同目录的 parameters.json 可核对原始材料、网格、频率列表、源表达式、孔隙元数据及文件来源。模型清单中的M编号与下面附录对应。原始含场解文件保留在原工作目录，本交付包中的MPH均应显示没有已计算解；这不是损坏，应执行“计算”生成自己的结果。')
paragraph('清理采用COMSOL的clearSolutionData接口并在重新打开文件后检查isEmpty，保留求解器序列。官方说明：https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_api_general.47.58.html。孔隙生成依据为本项目提供的张翔2011年论文及本地配置；详细来源与假设以本手册材料说明为准。')

heading('附录 本轮模型逐项索引')
paragraph('名称为 main_models 下的相对文件夹，文件名均为model.mph。频率按表中顺序对应freq_id的1、2、3；mode_id在全部主模型中均为0=S0、1=A0。hmax为最大网格控制尺寸，不是孔径。')
rows=[];csvrows=[]
for index,item in enumerate(main,1):
    source=Path(item['source']);record=read(source.parent/'results.json');sig=record['signature']
    name=Path(item['destination']).parent.as_posix().removeprefix('main_models/')
    count=record.get('finite_voids')
    meta=source.parent/'geometry/metadata.json'
    if count is None and meta.exists():count=read(meta)['pore_count']
    if count is None:count=0 if 'reference' in name else (1 if 'single' in name or 'maximum' in name else '见模型')
    frequencies='/'.join(f'{f/1000:g}' for f in sig['frequencies'])
    rows.append([f'M{index:03}',name,frequencies,count,f'{sig["mesh"]*1e6:g}'])
    csvrows.append(dict(编号=f'M{index:03}',相对路径=item['destination'],频率kHz=frequencies,孔数=count,最大网格微米=sig['mesh']*1e6))
table(['编号','模型文件夹','频率 kHz','孔数','hmax μm'],rows,[1.8,8.0,3.2,1.8,2.7])
with (DEST/'主模型参数索引.csv').open('w',encoding='utf-8-sig',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(csvrows[0]));writer.writeheader();writer.writerows(csvrows)
doc.core_properties.title='薄板孔隙群数值实验流程与复核手册'
doc.core_properties.subject='模型参数 操作流程 参考数据 验证状态'
doc.save(DOCX)
print(DOCX)
