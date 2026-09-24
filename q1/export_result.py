"""Package the current six-parameter fit without recalibration or data smoothing."""
from pathlib import Path
import csv, json, shutil, hashlib, html
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'q1'/'_result'
SRC = ROOT/'results'/'q1_v06'
if Path(__file__).resolve().parent.name=='01_model':
    OUT=Path(__file__).resolve().parent.parent
    ROOT=OUT/'01_model'/'repro_project'
    SRC=OUT/'06_raw_results'
VERIFY = SRC/'verification_phase'
for name in ['01_model','02_parameters','03_tables','04_figures/main','04_figures/variables','05_validation','06_raw_results','07_qa']:
    (OUT/name).mkdir(parents=True, exist_ok=True)

def read(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def write(path, rows):
    if not rows: return
    with Path(path).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def numeric(path):
    rows=read(path)
    return {k:np.array([float(r[k]) for r in rows]) for k in rows[0]}

def cp(src,dest):
    if Path(src).resolve()!=Path(dest).resolve():shutil.copy2(src,dest)

# Preserve all current six-parameter simulation states, grids and solver output.
if VERIFY.resolve()!=(OUT/'06_raw_results'/'verification_phase').resolve():
    shutil.copytree(VERIFY,OUT/'06_raw_results'/'verification_phase',dirs_exist_ok=True)
for name in ['fit_phase.mat','checkpoint_phase.mat','input_audit.csv']:
    if (SRC/name).exists(): cp(SRC/name,OUT/'06_raw_results'/name)
for name in ['grid_summary.csv','numerical_checks.csv','optimization_history.csv','phase_assumption_sensitivity.csv','verification_report.md']:
    cp(VERIFY/name,OUT/'05_validation'/name)
for name in ['q1_checks.csv']:
    if (ROOT/'results'/'checks'/name).exists():cp(ROOT/'results'/'checks'/name,OUT/'05_validation'/name)
if (ROOT/'results'/'q1_v05'/'local_phase_checks.csv').exists():
    cp(ROOT/'results'/'q1_v05'/'local_phase_checks.csv',OUT/'05_validation'/'local_phase_checks.csv')
for file in VERIFY.rglob('*.csv'):
    target=OUT/'03_tables'/'original'/file.relative_to(VERIFY)
    target.parent.mkdir(parents=True,exist_ok=True);cp(file,target)

fit=loadmat(SRC/'fit_phase.mat',simplify_cells=True)['fit']
cfg=json.loads((VERIFY/'fine'/'config.json').read_text(encoding='utf-8'))
cfg['calibrationCheckpoint']=''
(OUT/'02_parameters'/'config.json').write_text(json.dumps(cfg,ensure_ascii=False,indent=2),encoding='utf-8')
cp(VERIFY/'fitted_parameters.csv',OUT/'02_parameters'/'fitted_parameters.csv')
write(OUT/'02_parameters'/'all_config.csv',[{'name':k,'value':json.dumps(v,ensure_ascii=False)} for k,v in cfg.items()])
layernames=['Anode plate','Anode GDL','Anode CL','PEM','Cathode CL','Cathode GDL','Cathode plate']
write(OUT/'02_parameters'/'layer_parameters.csv',[{'layer':i+1,'name':n,'thickness_m':cfg['dxLayer'][i],
 'heat_capacity_J_m3_K':cfg['CLayer'][i],'dry_conductivity_W_m_K':cfg['kLayer'][i],
 'porosity':cfg['epsLayer'][i],'ionomer_fraction':1 if i==3 else (.3 if i in [2,4] else 0),
 'permeability_m2':cfg['K0Layer'][i],'contact_angle_deg':cfg['thetaLayer'][i]} for i,n in enumerate(layernames)])

# Self-contained MATLAB project; its root resolver remains valid after moving it.
project=OUT/'01_model'/'repro_project'; (project/'q1').mkdir(parents=True,exist_ok=True)
core=['config','project_root','read_data','make_grid','make_index','initial','unpack','properties','voltage',
      'face_fluxes','phase_rates','rhs','rhs_diagnostic','events','balance','observe','simulate',
      'fit_config','residual','calibrate','save_result','export_results','export_tables','verify_transport','plot_mesh']
for n in core: cp(ROOT/'q1'/f'q1_{n}.m',project/'q1'/f'q1_{n}.m')
for n in ['附件1.xlsx','附件2.xlsx']: cp(ROOT/n,project/n)
cp(SRC/'fit_phase.mat',project/'fit_phase.mat')
(project/'run_model.m').write_text("""% Fixed fitted parameters; coarse/medium/fine and time verification, no refit.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'q1'));
report=q1_verify_transport(root,'six_parameter',fullfile(root,'fit_phase.mat'),fullfile(root,'reproduced'));
""",encoding='utf-8')
(project/'recalibrate.m').write_text("""% Optional NEW six-parameter calibration; does not overwrite the delivered fit.
root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'q1'));
s=load(fullfile(root,'fit_phase.mat'),'fit');c=s.fit.cfg;
c.calibrationCheckpoint=fullfile(root,'new_checkpoint.mat');
[cases,~]=q1_read_data(root,c);
fit=q1_calibrate(cases,c,'coarse',300,1,1:6);
save(fullfile(root,'new_fit.mat'),'fit','cases');
""",encoding='utf-8')

data=[]
derived={
 'water_total_kgm2':('Total resident water','kg/m2','所有储水库总和，含初始束缚水'),
 'ice_total_kgm2':('Total frozen water','kg/m2','孔隙冰与膜内/离聚物冻结水总和'),
 'ice_total_mgcm2':('Total frozen water','mg/cm2','总冻结水面密度'),
 'pore_ice_mgcm2':('Pore ice','mg/cm2','孔隙冰面密度'),
 'ion_ice_mgcm2':('Ionomer / membrane ice','mg/cm2','膜内/离聚物冻结水面密度'),
 'resident_water_frozen_pct':('Frozen fraction of resident water','%','总冻结水/当前全部储水质量，不是局部体积分数'),
 'voltage_residual_V':('Voltage residual','V','模型减实验'),
 'temperature_residual_C':('Temperature residual','deg C','七层平均模型温度减实验'),
}
for k in [1,2]:
    d=numeric(VERIFY/'fine'/f'case{k}_observations.csv')
    d['water_total_kgm2']=sum(d[n] for n in ['water_bound_kgm2','water_vapor_kgm2','water_liquid_kgm2','water_pore_ice_kgm2','water_ion_ice_kgm2'])
    d['ice_total_kgm2']=d['water_pore_ice_kgm2']+d['water_ion_ice_kgm2']
    d['ice_total_mgcm2']=100*d['ice_total_kgm2']
    d['pore_ice_mgcm2']=100*d['water_pore_ice_kgm2'];d['ion_ice_mgcm2']=100*d['water_ion_ice_kgm2']
    d['resident_water_frozen_pct']=100*d['ice_total_kgm2']/d['water_total_kgm2']
    d['voltage_residual_V']=d['V_sim']-d['V_exp'];d['temperature_residual_C']=d['T_seven_C']-d['T_exp_C']
    assert len(d['time_s'])==184 and np.all(np.diff(d['time_s'])>0)
    assert np.allclose(d['ice_total_kgm2'],d['water_pore_ice_kgm2']+d['water_ion_ice_kgm2'])
    write(OUT/'03_tables'/f'case{k}_all_observables.csv',[{n:float(v[i]) for n,v in d.items()} for i in range(184)])
    data.append(d)

meta={
'V_sim':('Voltage','V','模拟电压'),'V_exp':('Observed voltage','V','实验电压'),
'T_seven_C':('Seven-layer mean temperature','deg C','七层厚度加权平均温度，标定口径'),
'T_MEA_C':('MEA mean temperature','deg C','五层MEA厚度加权平均温度'),
'T_exp_C':('Observed temperature','deg C','实验温度'),
'ice_volume_max':('Maximum local pore-ice volume fraction','1','空间最大 mi/rhoI，以局部单元总体积为分母'),
'ice_saturation_max':('Maximum pore-ice saturation','1','空间最大 mi/(rhoI*porosity)'),
'liquid_saturation_max':('Maximum liquid saturation','1','空间最大 ml/(rhoL*porosity)'),
'ice_peak_layer':('Layer of maximum pore ice','layer','无孔隙冰时为NaN'),
'ice_peak_x_um':('Position of maximum pore ice','um','无孔隙冰时为NaN'),
'lambda_aCL':('Anode CL hydration','1','阳极催化层平均含水量'),
'lambda_PEM':('Membrane hydration','1','膜平均含水量'),
'lambda_cCL':('Cathode CL hydration','1','阴极催化层平均含水量'),
'lambda_oversat_CL':('Maximum CL hydration excess','1','CL内最大lambda-lambdaSat'),
'CL_frozen_ion_water':('CL frozen ionomer water','kg/m2','催化层离聚物冻结水面密度'),
'gas_porosity_min':('Minimum gas-filled porosity','1','多孔层最小剩余气相孔隙率'),
'V_Erev':('Reversible potential','V','可逆电压'),'V_act':('Activation overpotential','V','活化损失'),
'V_ohm':('Ohmic overpotential','V','欧姆损失'),'V_con':('Concentration overpotential','V','浓差损失'),
'R_aCL':('Anode CL resistance','ohm m2','阳极催化层面积比电阻'),
'R_PEM':('Membrane resistance','ohm m2','膜面积比电阻'),
'R_cCL':('Cathode CL resistance','ohm m2','阴极催化层面积比电阻'),
'jlim_H':('Hydrogen limiting current density','A/m2','氢气极限电流密度'),
'jlim_O':('Oxygen limiting current density','A/m2','氧气极限电流密度'),
'reaction_weight_max':('Maximum reaction weight','1','阴极最大局部反应分配权重'),
'nucleation_progress':('Maximum nucleation progress','1','局部成核进度的空间最大值'),
'j_Acm2':('Applied current density','A/cm2','附件总电流除以303 cm2'),
'water_residual':('Water conservation residual','kg/m2','库存变化减生成水加外排水'),
'energy_residual':('Energy conservation residual','J/m2','焓变化减电化学热加对流及蒸气外排能量'),
}
for n in ['bound','vapor','liquid','pore_ice','ion_ice']:
    meta['water_'+n+'_kgm2']=(n.replace('_',' ').capitalize()+' water','kg/m2','按厚度积分的水库面密度')
meta.update(derived)
assert set(data[0])-{'time_s'}==set(meta)
write(OUT/'03_tables'/'variable_dictionary.csv',[{'field':'time_s','label':'Time','unit':'s','definition':'原始采样时间'}]+[
 {'field':k,'label':v[0],'unit':v[1],'definition':v[2]} for k,v in meta.items()])
write(OUT/'03_tables'/'final_water_inventory.csv',[{'case_C':-20-5*k,
 **{n:float(d[n][-1]) for n in ['water_bound_kgm2','water_liquid_kgm2','water_pore_ice_kgm2','water_ion_ice_kgm2','water_total_kgm2','ice_total_kgm2','resident_water_frozen_pct']}} for k,d in enumerate(data)])

plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':9,'axes.titlesize':11,'axes.labelsize':9,
 'legend.fontsize':8,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,'svg.fonttype':'none',
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'savefig.dpi':300})
colors=['#2878B5','#C96B32']; gallery=[]

def save(fig,name,title,description,sub='main'):
    path=OUT/'04_figures'/sub/name
    for ext in ['png','pdf','svg']:fig.savefig(path.with_suffix('.'+ext),facecolor='white',dpi=300)
    plt.close(fig)
    gallery.append({'name':name,'title':title,'description':description,'path':str(path.relative_to(OUT)).replace('\\','/')})

def decorate(ax,label):
    ax.set(xlabel='Time (s)',ylabel=label,xlim=(0,36.6)); ax.grid(alpha=.15)

def pair(title,height=3.5):
    fig,axes=plt.subplots(1,2,figsize=(9,height))
    fig.subplots_adjust(left=.09,right=.98,bottom=.25,top=.81,wspace=.33)
    fig.suptitle(title,y=.98,fontsize=12)
    for k,ax in enumerate(axes):ax.set_title(f'{-20-5*k} deg C')
    return fig,axes

# Figure contract: fit evidence, compartment interpretation, threshold dependency,
# and numerical verification are separate quantitative figures. All 184 samples
# are retained. No smoothing, inferred error bars or independent-validation claim.
for name,field,obs,title in [('01_voltage','V_sim','V_exp','Voltage fit'),('02_temperature','T_seven_C','T_exp_C','Temperature fit')]:
    fig,axes=pair(title)
    for ax,d in zip(axes,data):
        ax.plot(d['time_s'],d[obs],'.',color='#333333',ms=2,label='Observed')
        ax.plot(d['time_s'],d[field],color=colors[0],lw=1.3,label='Six-parameter fit')
        decorate(ax,meta[field][0]+' ('+meta[field][1]+')')
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=2,bbox_to_anchor=(.5,.015))
    save(fig,name,title,'两组工况各184点，均参与联合标定；无独立验证集。')

fig,axes=plt.subplots(1,3,figsize=(12,3.5));fig.subplots_adjust(left=.07,right=.98,bottom=.26,top=.80,wspace=.4)
fig.suptitle('Frozen-water compartments',y=.98,fontsize=12)
for ax,field in zip(axes,['pore_ice_mgcm2','ion_ice_mgcm2','ice_total_mgcm2']):
    for k,d in enumerate(data):ax.plot(d['time_s'],d[field],color=colors[k],ls=['-','--'][k],label=f'{-20-5*k} deg C')
    decorate(ax,'Frozen water (mg/cm2)');ax.set_title(meta[field][0]);ax.set_ylim(bottom=0,top=None)
    if all(np.all(d[field]==0) for d in data):
        ax.grid(False);ax.set_ylim(0,.001);ax.text(.5,.55,'No pore ice in this window',transform=ax.transAxes,ha='center',fontsize=8)
fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=2,bbox_to_anchor=(.5,.015))
save(fig,'03_ice_compartments','孔隙冰、膜内冻结水与总冻结水','总冰与膜内冰重合，因为本次基准方案孔隙冰为零。')

fig,axes=pair('Resident-water inventories',4)
poolcolors=['#2878B5','#C96B32','#5B8E67','#8A679E','#777777']
fields=['water_bound_kgm2','water_liquid_kgm2','water_ion_ice_kgm2','water_pore_ice_kgm2','water_vapor_kgm2']
for ax,d in zip(axes,data):
    for j,n in enumerate(fields):ax.plot(d['time_s'],100*d[n],label=n.removeprefix('water_').removesuffix('_kgm2').replace('_',' '),color=poolcolors[j],lw=1.2,ls=['-','-','-','--',':'][j])
    decorate(ax,'Water inventory (mg/cm2)');ax.set_ylim(bottom=0)
fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=5,bbox_to_anchor=(.5,.01))
save(fig,'04_water_inventories','水分配曲线','束缚水、液态水、膜内冻结水、孔隙冰和蒸气，含初始储水。')

fig,axes=pair('Liquid saturation and the nucleation threshold')
for ax,d in zip(axes,data):
    ax.plot(d['time_s'],100*d['liquid_saturation_max'],color=colors[0],label='Maximum liquid saturation')
    ax.axhline(3,color=colors[1],ls='--',label='Nucleation threshold: 3%')
    decorate(ax,'Liquid saturation (%)');ax.set_ylim(0,3.5)
fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=2,bbox_to_anchor=(.5,.015))
save(fig,'05_nucleation_threshold','液态水饱和度与成核阈值','两工况的最大液态水饱和度均低于3%，孔隙成核未启动。')

fig,axes=pair('Voltage decomposition',4)
for ax,d in zip(axes,data):
    for j,n in enumerate(['V_Erev','V_act','V_ohm','V_con']):ax.plot(d['time_s'],d[n],color=poolcolors[j],label=meta[n][0])
    decorate(ax,'Potential / overpotential (V)')
fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=2,bbox_to_anchor=(.5,0))
save(fig,'06_voltage_decomposition','电压损失分解','V=Erev-activation-ohmic-concentration。')

for field,label,suffix in [('V_sim','Voltage (V)','voltage'),('T_seven_C','Temperature (deg C)','temperature')]:
    fig,axes=pair('Fixed-parameter grid verification: '+suffix)
    for k,ax in enumerate(axes):
        for j,level in enumerate(['coarse','medium','fine']):
            d=numeric(VERIFY/level/f'case{k+1}_observations.csv')
            ax.plot(d['time_s'],d[field],color=poolcolors[j],ls=[':', '--','-'][j],label=level.capitalize())
        decorate(ax,label)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=3,bbox_to_anchor=(.5,.015))
    save(fig,'07_grid_'+suffix,'网格检查：'+suffix,'固定同一六参数，仅改变空间网格；粗31、中62、细124单元。')

st=read(VERIFY/'phase_assumption_sensitivity.csv')
fig,axes=plt.subplots(1,3,figsize=(12,3.8));fig.subplots_adjust(left=.07,right=.98,bottom=.30,top=.81,wspace=.38)
fig.suptitle('Phase-assumption sensitivity (fixed fitted parameters)',y=.98,fontsize=12)
for ax,field,label in zip(axes,['V_RMSE','T_RMSE_C','ice_peak'],['Voltage RMSE (V)','Temperature RMSE (deg C)','Peak pore-ice volume fraction']):
    for k in [0,1]:
        r=st[k::2];ax.plot(np.arange(5),[float(x[field]) for x in r],marker='o',ms=4,color=colors[k],ls=['-','--'][k],label=f'{-20-5*k} deg C')
    ax.set_xticks(np.arange(5),['.02 / 1','.03 / .5','.03 / 1','.03 / 2','.04 / 1'],rotation=25,ha='right',rotation_mode='anchor')
    ax.set(xlabel='Threshold / freezing rate (1/s)',ylabel=label);ax.grid(alpha=.15)
fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=2,bbox_to_anchor=(.5,.005))
save(fig,'08_phase_sensitivity','相变假设敏感性','中网格固定已标定参数；成核阈值改变冰量，V/T拟合误差变化很小。')

hist=read(VERIFY/'optimization_history.csv')
fig,ax=plt.subplots(figsize=(7,3.5));fig.subplots_adjust(left=.13,right=.97,bottom=.18,top=.84)
ax.plot([int(r['funccount']) for r in hist],[float(r['weighted_SSE']) for r in hist],'-o',ms=3,color=colors[0])
ax.set(xlabel='Function evaluations',ylabel='Calibration weighted SSE',title='Six-parameter optimization history');ax.grid(alpha=.15)
save(fig,'09_optimization','标定迭代过程','残差尺度为0.02 V和0.1 ℃；最后因步长容差停止，不代表全局最优。')

# Exhaustive scalar output atlas: every original observable plus derived water totals.
for field,(label,unit,definition) in meta.items():
    fig,ax=plt.subplots(figsize=(7,3.7));fig.subplots_adjust(left=.15,right=.97,bottom=.26,top=.83)
    ax.set_title(label)
    anyvalid=False
    for k,d in enumerate(data):
        valid=np.isfinite(d[field]);anyvalid|=bool(valid.any())
        ax.plot(d['time_s'],d[field],color=colors[k],ls=['-','--'][k],lw=1.3,label=f'{-20-5*k} deg C')
    decorate(ax,unit)
    if not anyvalid:
        ax.grid(False);ax.set_yticks([])
        ax.text(.5,.5,'Undefined: no pore ice in this window',ha='center',transform=ax.transAxes,fontsize=9)
    if anyvalid and all(np.all(d[field]==0) for d in data):
        ax.grid(False)
        ax.set_ylim(0,1);ax.text(.5,.5,'Identically zero over the observation window',ha='center',transform=ax.transAxes,fontsize=8)
    fig.legend(*ax.get_legend_handles_labels(),loc='lower center',ncol=2,bbox_to_anchor=(.5,.015))
    save(fig,field,label,definition,'variables')

write(OUT/'04_figures'/'figure_index.csv',gallery)
cp(Path(__file__),OUT/'01_model'/'export_result.py')
contract='''# Figure and data contract\n\nOnly the current six-parameter solution is included.\nAll 184 original times per condition are retained; no smoothing, exclusions, simulated error bars or new optimization.\nFigures 01–02 establish the training fit, 03–05 explain compartment/threshold behavior, 07 checks numerical precision, and 08 tests the unmeasured phase assumptions. The variable atlas documents every output, without treating each diagnostic as independent scientific evidence.\nPython/matplotlib; 7–12 inch width; text 8–12 pt; white background; editable PDF/SVG and 200 dpi PNG.\nBoth conditions are training data. Ice has no experimental validation.\nLegacy MATLAB graphics are retained unchanged under raw results; newly generated graphics are in 04_figures.\n'''
(OUT/'07_qa'/'figure_contract.md').write_text(contract.replace('200 dpi','300 dpi'),encoding='utf-8')
cards=''.join(f'<article><h2>{html.escape(g["title"])}</h2><p>{html.escape(g["description"])}</p><a href="{g["path"]}.pdf">PDF</a> · <a href="{g["path"]}.svg">SVG</a><br><img loading="lazy" src="{g["path"]}.png"></article>' for g in gallery)
(OUT/'index.html').write_text('''<!doctype html><html lang="zh"><meta charset="utf-8"><title>Q1 六参数结果</title><style>body{font:16px system-ui;max-width:1100px;margin:36px auto;padding:0 24px;background:#f5f6f8;color:#172b3a}article{background:white;padding:24px;margin:24px 0;border-radius:10px}img{max-width:100%}a{color:#2878b5}h1{font-size:32px}h2{font-size:21px}</style><h1>Q1 当前六参数标定结果</h1><p>两组实验均用于标定；孔隙冰为零不等于总冻结水为零。</p><p><a href="README.md">文件导航</a> · <a href="01_model/完整模型说明.md">完整模型说明</a> · <a href="02_parameters/标定结果.md">参数与指标</a> · <a href="03_tables/variable_dictionary.csv">变量字典</a></p>'''+cards+'</html>',encoding='utf-8')
print(json.dumps({'output':str(OUT),'figures':len(gallery),'samples_per_case':184,'variables':len(meta)},ensure_ascii=False))
