function result=q1_compare_revision(root,out)
% Keep original and revised residual scales; never optimize reference ice.
base=fullfile(root,'results','q1_v06');rows={};
modes={'direct','transport','phase'};
for z=1:3
    mode=modes{z};npar=5+(z==3);
    if z==1
        t=out.baseline;
        [cases,~]=q1_read_data(root);
        old=zeros(2,1);new=old;
        for k=1:2
            n=numel(cases(k).t);
            old(k)=n*((t.V_RMSE(k)/.02)^2+(t.T_RMSE_C(k)/.5)^2);
            new(k)=n*((t.V_RMSE(k)/.02)^2+(t.T_RMSE_C(k)/.1)^2);
        end
        mesh=NaN;time=NaN;flag=NaN;
    else
        r=out.(mode);t=readtable(fullfile(base,['verification_' mode],'fine','summary.csv'));
        f=r.summary(r.summary.grid=="fine",:);
        old=f.weighted_SSE;new=f.calibration_SSE;
        mesh=r.meshPass;time=r.timePass;flag=r.fit.exitflag;
    end
    for k=1:height(t)
        rows(end+1,:)={string(mode),npar,string(t.case_name(k)),t.V_RMSE(k), ...
            t.T_RMSE_C(k),old(k),new(k),t.ice_onset_s(k),t.ice_volume_peak(k), ...
            t.water_error(k),t.energy_error(k),mesh,time,flag}; %#ok<AGROW>
    end
end
result.table=cell2table(rows,'VariableNames',{'model','n_parameters','case_name', ...
    'V_RMSE','T_RMSE_C','weighted_SSE','calibration_SSE','ice_onset_s', ...
    'ice_volume_peak','water_error','energy_error','mesh_pass','time_pass','exitflag'});
writetable(result.table,fullfile(base,'model_comparison.csv'));
t=result.table;a=t(t.model=="transport",:);b=t(t.model=="phase",:);
result.relativeImprovement=1-sum(b.calibration_SSE)/sum(a.calibration_SSE);
% A transparent working rule, not a statistical significance test.
eligible=[out.transport.meshPass && out.transport.timePass && out.transport.fit.exitflag>0, ...
    out.phase.meshPass && out.phase.timePass && out.phase.fit.exitflag>0];
result.selectedMode='transport';
if eligible(2) && (~eligible(1) || result.relativeImprovement>.05)
    result.selectedMode='phase';
elseif ~any(eligible) && result.relativeImprovement>.05
    result.selectedMode='phase';
end
result.finalSelection=any(eligible);
result.selectionScope='voltage_temperature_calibration';
fid=fopen(fullfile(base,'修订结果说明.md'),'w','n','UTF-8');
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# Q1 五参数与六参数比较\n\n');
fprintf(fid,'使用当前目录附件2的两组电压与温度联合标定，共用一套参数；已知几何和热物性保持原值。以下均为训练误差。\n\n');
fprintf(fid,'## 模型与参考图对应\n\n');
fprintf(fid,'一维有限体积水热电化学耦合；五层 MEA 外加两侧板材，共七层计算域。保留界面通量连续、三相水迁移、局部成核、相变潜热、冰占孔和反应面积反馈。采用 ode15s 隐式积分，催化层和膜区域分层加密。\n\n');
fprintf(fid,'温度拟合口径沿用七层厚度加权平均 T_seven_C；五层 MEA 平均 T_MEA_C 同时输出，但未替换为标定温度。\n\n');
fprintf(fid,'参考图只提供建模流程，没有可核实的数值目标。因此没有将先前代码中硬编码的参考 RMSE、起冰时间和冰量作为标定数据或已核实比较值。\n\n');
if isfile(fullfile(base,'reference_model.png'))
    fprintf(fid,'用户提供的参考流程图保存在 [reference_model.png](reference_model.png)。\n\n');
end
fprintf(fid,'## 同一细网格上的比较\n\n');
fprintf(fid,'| 模型 | 工况 | 电压 RMSE/V | 温度 RMSE/℃ | 旧权重 SSE | 新权重 SSE | 起冰/s | 峰值冰体积分数 |\n|---|---|---:|---:|---:|---:|---:|---:|\n');
for k=1:height(t)
    onset="未达阈值";
    if isfinite(t.ice_onset_s(k)),onset=string(sprintf('%.1f',t.ice_onset_s(k)));end
    fprintf(fid,'| %s | %s | %.6f | %.6f | %.3f | %.3f | %s | %.7f |\n', ...
        t.model(k),t.case_name(k),t.V_RMSE(k),t.T_RMSE_C(k),t.weighted_SSE(k), ...
        t.calibration_SSE(k),onset,t.ice_volume_peak(k));
end
fprintf(fid,'\n旧权重残差尺度为 0.02 V/0.5 ℃；新权重为 0.02 V/0.1 ℃。这些是拟合权衡，不是测量不确定度。禁止跨权重直接比较目标值。\n\n');
fprintf(fid,'六参数相对五参数的细网格新权重 SSE 改善为 %.2f%%。\n\n',100*result.relativeImprovement);
for k=1:height(b)
    if isnan(b.ice_onset_s(k))
        fprintf(fid,'六参数 %s 工况在 0–36.6 s 内未达到孔隙起冰判据（峰值 %.8g）。因此该工况的电压/温度拟合改善不能用来验证孔隙冰堵或冻结速率。\n\n', ...
            b.case_name(k),b.ice_volume_peak(k));
    end
end
fprintf(fid,'工作选择规则：优先通过空间与时间精度检查且优化正常停止的方案；两者均合格时，六参数须降低至少 5%% 的新权重 SSE 才抵偿增加参数的复杂度。这是预设工作阈值，不是显著性检验。\n\n');
if result.finalSelection
    fprintf(fid,'按上述规则选择 `%s` 作为电压/温度拟合方案；该选择不代表结冰机制已通过实验验证。\n\n',result.selectedMode);
else
    fprintf(fid,'两种方案各自均未满足全部优化停止与数值条件，不能宣布最终模型；暂以 `%s` 为后续检查候选。\n\n',result.selectedMode);
end
for mode={'transport','phase'}
    r=out.(mode{1});
    fprintf(fid,'- %s：退出标志 %d；中→细网格通过 %d；细网格时间精度通过 %d；粗网格优化目标 %.6f。\n', ...
        mode{1},r.fit.exitflag,r.meshPass,r.timePass,r.fit.resnorm);
    medium=r.summary(r.summary.grid=="medium",:);
    fine=r.summary(r.summary.grid=="fine",:);
    for k=1:height(fine)
        fprintf(fid,'  %s 的中/细网格冰峰值：%.8g / %.8g。\n', ...
            fine.case_name(k),medium.ice_volume_peak(k),fine.ice_volume_peak(k));
    end
end
fprintf(fid,'\n数值阈值为最大电压差 5 mV、温度差 0.05 ℃、冰体积分数差 0.005；守恒另按水 1e-4、能量 1e-3 检查。若网格未通过，需加密并视情况重新标定。优化正常停止不代表全局最优。\n\n');
fprintf(fid,'冰量检查使用绝对容差；当冰峰值很小时，通过该检查并不代表冰量的相对误差或起冰时刻已经收敛。应结合上列峰值变化和敏感性结果解释。\n\n');
fprintf(fid,'## 结冰假设与证据边界\n\n');
fprintf(fid,'成核阈值 0.03、冻结速率 1 s⁻¹ 为假设。固定已标定参数，分别试算阈值 0.02/0.03/0.04 和冻结速率 0.5/1/2 s⁻¹，见候选方案 verification_* 内 phase_assumption_sensitivity.csv。\n\n');
fprintf(fid,'起冰为每 0.2 s 输出的多孔层孔隙冰空间最大 mi/rhoI 首次超过 1e-6；峰值为时空最大 mi/rhoI，并非孔隙冰饱和度，也不包含离聚物/膜内冻结水。后者由 water_ion_ice_kgm2 单独输出，膜内仍采用自身含水容量冻结律。冰量无实验观测，不能据此宣称冰量预测已验证。\n\n');
fprintf(fid,'两种温度均参与拟合，不构成独立验证。参考图要求的独立预测检验仍需留出工况或新实验。接近边界的参数见各 verification_* 内 fitted_parameters.csv，不应解释为已独立识别的材料常数。\n');
selected=out.(result.selectedMode);
if isfield(selected,'sensitivity')
    st=selected.sensitivity;
    fprintf(fid,'\n## 相变假设敏感性结果（中网格）\n\n');
    fprintf(fid,'全部标定参数固定，仅改变表列假设，未重新拟合。\n\n');
    fprintf(fid,'| 阈值 | 冻结速率/s⁻¹ | 工况 | 电压 RMSE/V | 温度 RMSE/℃ | 起冰/s | 冰峰值 |\n|---:|---:|---|---:|---:|---:|---:|\n');
    for k=1:height(st)
        onset="未达阈值";
        if isfinite(st.ice_onset_s(k)),onset=string(sprintf('%.1f',st.ice_onset_s(k)));end
        fprintf(fid,'| %.2f | %.1f | %s | %.6f | %.6f | %s | %.8g |\n', ...
            st.nucleation_threshold(k),st.kfreeze_s_inv(k),st.case_name(k), ...
            st.V_RMSE(k),st.T_RMSE_C(k),onset,st.ice_peak(k));
    end
    m=abs(st.nucleation_threshold-.03)<1e-12;
    if any(m) && all(st.ice_peak(m)==0)
        fprintf(fid,'\n基准阈值 0.03 下，冻结速率 0.5/1/2 s⁻¹ 均未产生孔隙冰，这组电压/温度观测不能辨识该冻结速率。成核阈值改变时是否起冰，应作为模型假设依赖报告，不能解释为已有冰量实验支持。\n');
    end
end
disp(result.table);
end
