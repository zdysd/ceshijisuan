"""Chinese layout diagram using actual pores from the 48 mm thin-plate case."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Polygon

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'output/plate_layout';OUT.mkdir(exist_ok=True)
DATA=ROOT/'output/guided_collective'
s=json.loads((DATA/'study_settings.json').read_text(encoding='utf-8'))
g=json.loads((DATA/'d80/L48/seed_101/geometry/geometry.json').read_text(encoding='utf-8'))
start=s['start']*1000;end=start+48;sponge=s['sponge']*1000;total=s['length']*1000
blue='#d9eafa';yellow='#fff0c9';edge='#37506b';teal='#087f8c';red='#b5542a'
plt.rcParams.update({'font.family':'Microsoft YaHei','axes.unicode_minus':False,'font.size':11,'svg.fonttype':'none'})
fig=plt.figure(figsize=(15,7.7),facecolor='white')
grid=fig.add_gridspec(2,2,height_ratios=[2.3,1.15],hspace=.34,wspace=.18)
ax=fig.add_subplot(grid[0,:]);left=fig.add_subplot(grid[1,0]);right=fig.add_subplot(grid[1,1])
fig.suptitle('二维薄板中无孔区与孔隙群的布置',fontsize=20,fontweight='bold',y=.98)
fig.text(.5,.925,'示例：厚度 1 毫米  ·  孔隙段长 48 毫米  ·  80 微米孔径组  ·  随机种子 101',ha='center',fontsize=12,color='#435366')


def draw_board(a):
    a.add_patch(Rectangle((0,0),total,1,facecolor=blue,edgecolor='none'))
    a.add_patch(Rectangle((start,0),end-start,1,facecolor=yellow,edgecolor='none'))
    for j in range(70):
        x=sponge+(total-sponge)*j/70
        a.add_patch(Rectangle((x,0),(total-sponge)/70+.003,1,facecolor=plt.cm.Greys(.10+.31*j/69),edgecolor='none'))
    for pore in g['polygons']:
        p=np.array(pore['exterior'])/1000+[start,0]
        a.add_patch(Polygon(p,closed=True,facecolor='white',edgecolor='#302e2a',lw=.45))
        for ring in pore.get('interiors',[]):
            a.add_patch(Polygon(np.array(ring)/1000+[start,0],closed=True,facecolor=yellow,edgecolor='#302e2a',lw=.45))
    a.add_patch(Rectangle((0,0),total,1,fill=False,edgecolor=edge,lw=1.6))


draw_board(ax);ax.set_xlim(-6,total+18);ax.set_ylim(-1.03,3.1);ax.axis('off')
ax.text(0,2.78,'整体布局  厚度方向放大显示',fontweight='bold',fontsize=13)
ax.annotate('',xy=(66,2.12),xytext=(5,2.12),arrowprops=dict(arrowstyle='->',color=edge,lw=2))
ax.text(35,2.27,'波沿板长方向传播',ha='center',color=edge)
ax.plot([0,0],[0,1],color='#235fa5',lw=5)
ax.annotate('左端激励\nS0 或 A0',xy=(0,.55),xytext=(12,1.23),ha='center',fontsize=11,
            arrowprops=dict(arrowstyle='-',color='#235fa5'),color='#235fa5')
ax.text(42,.43,'无孔传播区',ha='center',fontsize=13,color=edge)
ax.text((start+end)/2,.45,'孔隙群',ha='center',fontsize=13,color='#79501e',bbox=dict(facecolor=yellow,edgecolor='none',pad=1))
ax.text((sponge+total)/2,.47,'吸收区\n逐渐衰减传播波',ha='center',va='center',fontsize=12,color='#434343')
ax.annotate('孔后无孔接收区',xy=((end+sponge)/2,.8),xytext=(190,1.32),ha='center',fontsize=11,
            arrowprops=dict(arrowstyle='-',color=edge),color=edge)
for x in np.array(s['pre_positions'])*1000:
    ax.plot([x,x],[-.08,1.08],ls='--',color=teal,lw=1)
post=np.array([end+8,end+10,end+12])
for x in post:ax.plot([x,x],[-.08,1.08],ls='--',color=red,lw=1)
ax.annotate('孔前接收截面\n约 80、82、84 毫米',xy=(82,1.08),xytext=(76,1.86),ha='center',fontsize=10,color=teal,
            arrowprops=dict(arrowstyle='-',color=teal))
ax.annotate('孔后接收截面\n约 148、150、152 毫米',xy=(150,1.08),xytext=(152,1.86),ha='center',fontsize=10,color=red,
            arrowprops=dict(arrowstyle='-',color=red))
for x in [start,end,sponge]:ax.plot([x,x],[-.18,1.1],color=edge,lw=.8,alpha=.7)
for x0,x1,label in [(0,start,'约 92 毫米'),(start,end,'48 毫米'),(end,sponge,'14 毫米'),(sponge,total,'约 184 毫米')]:
    ax.annotate('',xy=(x1,-.36),xytext=(x0,-.36),arrowprops=dict(arrowstyle='<->',color=edge,lw=1))
    ax.text((x0+x1)/2,-.50,label,ha='center',va='top',fontsize=10 if x1-x0>20 else 9,color=edge)
for x,label in [(0,'0'),(start,'92'),(end,'140'),(sponge,'154'),(total,'338')]:
    ax.text(x,-.85,label,ha='center',fontsize=10,color='#536170')
ax.text(total+6,-.85,'位置\n（毫米）',va='center',fontsize=9,color='#536170')
ax.annotate('',xy=(total+4,1),xytext=(total+4,0),arrowprops=dict(arrowstyle='<->',color=edge))
ax.text(total+6,.5,'厚度\n1 毫米',va='center',fontsize=10,color=edge)

for a,bounds,boundary,title in [
    (left,(start-2,start+6),start,'孔隙段起点局部放大'),
    (right,(end-3,end+5),end,'孔隙段终点局部放大')]:
    draw_board(a);a.set_xlim(*bounds);a.set_ylim(-.06,1.06);a.set_aspect('equal',adjustable='box')
    a.axvline(boundary,color='#946a2f',ls='--',lw=1)
    a.set_title(title+'  孔形与厚度按同一比例显示',fontsize=11,pad=14)
    a.set_xlabel('传播方向位置（毫米）',fontsize=10,labelpad=6)
    a.set_ylabel('厚度方向（毫米）',fontsize=10)
    a.set_yticks([0,.5,1]);a.tick_params(labelsize=9)
    a.spines[['top','right']].set_visible(False)
    if a is left:
        a.text(start-1,.5,'无孔区',ha='center',va='center',color=edge,fontsize=10)
        a.text(start+3,1.12,'含小孔的连续基体',ha='center',color='#79501e',fontsize=10)
    else:
        a.text(end+2.5,.5,'孔后无孔区',ha='center',va='center',color=edge,fontsize=10)
        a.text(end-1.5,1.12,'孔隙群',ha='center',color='#79501e',fontsize=10)

fig.text(.06,.074,'颜色说明：浅蓝为无孔基体；浅黄仅标出孔隙分布范围，基体材料相同；白色黑边轮廓为实际小孔。',fontsize=11,color='#435366')
fig.text(.06,.036,'孔隙段内共 141 个不规则孔，面积孔隙率为 1.5%；本图为传播方向与厚度方向的截面，尺寸标注作了四舍五入。',fontsize=10,color='#657180')
fig.subplots_adjust(top=.87,bottom=.18,left=.06,right=.965)
fig.savefig(OUT/'二维薄板孔隙段示意图.png',dpi=200,bbox_inches='tight',facecolor='white')
fig.savefig(OUT/'二维薄板孔隙段示意图.svg',bbox_inches='tight',facecolor='white')
print(OUT/'二维薄板孔隙段示意图.png')
