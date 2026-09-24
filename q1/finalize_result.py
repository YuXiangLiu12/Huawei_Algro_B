"""Write six-parameter documentation and verify exported numeric tables."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
from scipy.io import loadmat

out=Path(__file__).resolve().parent/'_result'
def rows(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
cfg=json.loads((out/'02_parameters/config.json').read_text(encoding='utf-8'))
fit=loadmat(out/'06_raw_results/fit_phase.mat',simplify_cells=True)['fit']
pars=rows(out/'02_parameters/fitted_parameters.csv')
metrics=rows(out/'06_raw_results/verification_phase/fine/summary.csv')
grid=rows(out/'05_validation/grid_summary.csv')
fine=[r for r in grid if r['grid']=='fine']
note=['# 当前六参数标定结果\n',
'本结果包仅包含当前六参数局部相变方案。几何和已知热物性保持原值。电压、温度两组实验均参与训练拟合；无冰量实验、无独立留出验证集。\n',
'## 标定参数\n','| 参数 | 数值 | 单位 | 下界 | 上界 | 触及边界 |','|---|---:|---|---:|---:|---|']
units=['A/m²','J/mol','1','1','1','s⁻¹']
names=['−20 ℃参考交换电流密度 j0,253','活化能 Ea','CL电导率修正 fCL','动力学含水量指数 β','束缚水扩散倍率 aD','束缚水释放速率 kbl']
for r,n,u in zip(pars,names,units):
    note.append(f"| {n} | {float(r['values']):.10g} | {u} | {float(r['lower']):g} | {float(r['upper']):g} | {'是' if r['at_bound']=='1' else '否'} |")
note += [f"\n298.15 K参考交换电流密度 j0ref={cfg['j0ref']:.10g} A/m²，是由j0,253与Ea换算的派生量，不是第七个独立参数。\n",
'固定成核阈值0.03、成核速率1、冻结速率1 s⁻¹；其它配置见 all_config.csv 和 config.json。活化能位于下界，不能解释为独立测得的材料常数。\n',
'## 细网格误差\n','| 工况 | 电压RMSE/V | 温度RMSE/℃ | 旧权重SSE | 实际标定权重SSE |','|---|---:|---:|---:|---:|']
for r in fine:note.append(f"| {r['case_name']} | {float(r['V_RMSE']):.8f} | {float(r['T_RMSE_C']):.8f} | {float(r['weighted_SSE']):.6f} | {float(r['calibration_SSE']):.6f} |")
note += [f"\n细网格新权重总SSE：{sum(float(r['calibration_SSE']) for r in fine):.9f}；旧权重总SSE：{sum(float(r['weighted_SSE']) for r in fine):.9f}。",
f"粗网格优化目标：{fit['resnorm']:.9f}；退出标志{int(fit['exitflag'])}，步长达到设定容差。不同网格或不同权重下的目标值必须分别标记。\n",
'新权重残差尺度为0.02 V/0.1 ℃，旧权重为0.02 V/0.5 ℃；不是测量不确定度。\n',
'## 冻结水的正确解释\n',
'细网格0–36.6 s内，两工况的孔隙冰体积分数均为0，孔隙起冰时间为未达到判据；膜内冻结水非零。\n',
'| 工况 | 末时刻孔隙冰/(mg/cm²) | 末时刻膜内冻结水/(mg/cm²) | 总冻结水/(mg/cm²) | 总冻结水/全部库存水 |','|---|---:|---:|---:|---:|']
audit=[]
for k in [1,2]:
    rr=rows(out/f'03_tables/case{k}_all_observables.csv');d=rr[-1]
    note.append(f"| {-15-5*k} ℃ | {float(d['pore_ice_mgcm2']):.8f} | {float(d['ion_ice_mgcm2']):.8f} | {float(d['ice_total_mgcm2']):.8f} | {float(d['resident_water_frozen_pct']):.6f}% |")
    assert len(rr)==184
    raw=rows(out/f'06_raw_results/verification_phase/fine/case{k}_observations.csv')
    for col in raw[0]:
        a=np.array([float(r[col]) for r in raw]);b=np.array([float(r[col]) for r in rr])
        assert np.array_equal(a,b,equal_nan=True),col
    vr=np.sqrt(np.mean([float(r['voltage_residual_V'])**2 for r in rr]))
    tr=np.sqrt(np.mean([float(r['temperature_residual_C'])**2 for r in rr]))
    assert np.isclose(vr,float(metrics[k-1]['V_RMSE']),rtol=1e-10)
    assert np.isclose(tr,float(metrics[k-1]['T_RMSE_C']),rtol=1e-10)
    audit.append({'case':k,'samples':184,'original_columns_unchanged':True,'rmse_recomputed':True})
note += ['\n这里的“总冻结水/全部库存水”包含初始束缚水，不能与“最大局部孔隙冰体积分数”互换。\n',
'## 数值验证与敏感性\n',
'粗、中、细网格与细网格时间加严均有效，预设网格/时间差异和水/能量守恒检查通过。详情在05_validation。成核阈值改为0.02时，−20 ℃/−25 ℃在中网格约29.6/27.4 s出现孔隙冰，峰值约0.000664/0.001816，而电压和温度拟合几乎不变。因此冰量和成核阈值并未被现有电压/温度数据确定。\n',
'原始table1/table2为0、5、…、35 s输出。T_error_pct使用摄氏温度绝对值作分母，仅保留为原程序输出，不应作为具有温标不变性的温度相对误差。解释拟合精度时优先使用温度绝对误差或RMSE。\n']
(out/'02_parameters/标定结果.md').write_text('\n'.join(note),encoding='utf-8')
readme='''# Q1 六参数结果包

本目录只保留当前六参数方案，不包含其它参数数量的结果。可以直接双击 **index.html** 浏览所有新生成的图。

## 文件导航

| 目录/文件 | 内容 |
|---|---|
| 01_model/完整模型说明.md | 计算域、完整控制方程、物性、相变与成核、电压、边界与初值、标定及数值方法、指标定义 |
| 01_model/repro_project | 自包含MATLAB源码、原始附件和六参数MAT文件；run_model.m验证，recalibrate.m可选重新标定 |
| 01_model/export_result.py | 用Python从本包原始结果重新生成全部图表；需numpy、scipy、matplotlib，不运行优化 |
| 02_parameters/标定结果.md | 六参数数值、边界、误差与冻结水解释 |
| 02_parameters/*.csv、config.json | 所有配置、逐层参数和六个拟合参数 |
| 03_tables/case1_all_observables.csv | −20 ℃、184时刻的所有原始观测量及派生总水/总冰指标 |
| 03_tables/case2_all_observables.csv | −25 ℃对应数据 |
| 03_tables/variable_dictionary.csv | 逐列定义与单位；NaN含义有说明 |
| 03_tables/original | 原程序所有CSV表，按粗中细网格和验证结果保留 |
| 04_figures/main | 10张主题汇总图：拟合、冻结水分区、水分配、成核、电压分解、网格、敏感性、优化历史 |
| 04_figures/variables | 43个标量输出/派生量的完整时间曲线，均含两种温度 |
| 04_figures/figure_index.csv | 图名与科学含义；每图提供PNG、PDF、SVG |
| 05_validation | 网格、时间、守恒、机制、优化历史与相变敏感性检查 |
| 06_raw_results | 本次六参数完整MAT状态、拟合与检查点、原始图表；原图保持不变 |
| 07_qa | 数据一致性、绘图检查、可视检查预览及文件清单 |

## 使用和口径

case1为−20 ℃，case2为−25 ℃。主交付轨迹采用细网格；敏感性表采用中网格，均清楚区分。
53张新图均保留全部184个原始时间点，无平滑、无删点、无虚构误差棒；个别位置/层号字段在无孔隙冰时为NaN，图中明确标为未定义。

孔隙冰曲线为0不代表总冰为0。查看03_ice_compartments及总冻结水派生列，可看到膜内冻结水。两种温度都参加了标定，因此拟合误差不是独立预测验证。

MAT文件保留完整空间—时间状态与求解器信息，文件夹体积较大。若只阅读，可从本README、模型说明、参数结果和index.html开始。
'''
(out/'README.md').write_text(readme,encoding='utf-8')
(out/'07_qa/data_checks.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
print('Documentation and all original data columns verified.')
