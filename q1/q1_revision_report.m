function report=q1_revision_report(root,mode)
% Compare recalibrated local phase model against the direct model and references.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
if nargin<2,mode='phase';end
base=fullfile(root,'results','q1_v06');
fitfile=fullfile(base,['fit_' mode '.mat']);folder=fullfile(base,['verification_' mode]);
report=q1_verify_transport(root,['local_' mode],fitfile,folder);
new=readtable(fullfile(folder,'fine','summary.csv'));
old=readtable(fullfile(root,'results','q1_v05','verification_direct','fine','summary.csv'));
refs=[.02787 .25182 NaN NaN;.03744 .10205 30.4 .03128];
rows={};
for k=1:2
    fields={'V_RMSE','T_RMSE_C','ice_onset_s','ice_volume_peak'};
    for j=1:numel(fields)
        rows(end+1,:)={string(new.case_name(k)),string(fields{j}), ...
            old.(fields{j})(k),new.(fields{j})(k),refs(k,j)}; %#ok<AGROW>
    end
end
comparison=cell2table(rows,'VariableNames',{'case_name','metric','previous_direct', ...
    'revised_local','reference_model'});
writetable(comparison,fullfile(folder,'model_comparison.csv'));disp(comparison);
% Nucleation threshold is an unmeasured assumption. Never fit to reference ice.
s=load(fitfile,'fit','cases');rows={};
settings=[.02 1;.03 .5;.03 1;.03 2;.04 1];
for z=1:size(settings,1)
    c=s.fit.cfg;c.nucleationThreshold=settings(z,1);c.kfreeze=settings(z,2);
    for k=1:2
        r=q1_simulate(s.cases(k),c,'medium');
        assert(strcmp(r.status,'VALID'),'Sensitivity run failed: %s',r.event);
        o=r.obs;ice=[o.ice_volume_max];idx=find(ice>1e-6,1);onset=NaN;
        if ~isempty(idx),onset=o(idx).time_s;end
        rows(end+1,:)={c.nucleationThreshold,c.kfreeze,string(s.cases(k).name), ...
            sqrt(mean(([o.V_sim]'-s.cases(k).V_exp).^2)), ...
            sqrt(mean(([o.T_seven_C]'-s.cases(k).T_exp_C).^2)),onset,max(ice)}; %#ok<AGROW>
    end
end
sensitivity=cell2table(rows,'VariableNames',{'nucleation_threshold','kfreeze_s_inv', ...
    'case_name','V_RMSE','T_RMSE_C','ice_onset_s','ice_peak'});
writetable(sensitivity,fullfile(folder,'phase_assumption_sensitivity.csv'));
report.comparison=comparison;report.sensitivity=sensitivity;
save(fullfile(folder,'revision_report.mat'),'report');
writeSummary(folder,report,new,old,refs,mode);
end

function writeSummary(folder,report,new,old,refs,mode)
fid=fopen(fullfile(folder,'修订结果说明.md'),'w','n','UTF-8');
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# Q1 局部液态水—冰相变修订\n\n');
fprintf(fid,'修订分支：%s。电压和温度采用附件原始数据进行联合标定；截图中的结冰指标来自其他模型，只作对比。\n\n',mode);
fprintf(fid,'## 模型改动\n\n');
fprintf(fid,'- 启用过冷液态水、液态水迁移和局部成核；过饱和束缚水先进入液态水，再由局部成核进度控制冻结。\n');
fprintf(fid,'- 保留相变潜热、局部冰占孔、有效扩散和反应面积反馈，保留五层 MEA 与两侧板材的七层计算域。\n');
fprintf(fid,'- 已知几何、材料热物性与换热系数保持原值；调整电化学、输运及束缚水释放速率参数。\n');
fprintf(fid,'- 温度残差尺度由 0.5 ℃ 改为 0.1 ℃，电压尺度保持 0.02 V。这是拟合权衡，不是测量不确定度。\n\n');
fprintf(fid,'## 细网格指标对比\n\n');
fprintf(fid,'| 工况 | 指标 | 原直接结冰模型 | 本次修订 | 截图参考模型 |\n|---|---|---:|---:|---:|\n');
for k=1:2
    fprintf(fid,'| %s | 电压 RMSE / V | %.6f | %.6f | %.5f |\n',string(new.case_name(k)),old.V_RMSE(k),new.V_RMSE(k),refs(k,1));
    fprintf(fid,'| %s | 温度 RMSE / ℃ | %.6f | %.6f | %.5f |\n',string(new.case_name(k)),old.T_RMSE_C(k),new.T_RMSE_C(k),refs(k,2));
end
fprintf(fid,'| −25 ℃ | 起冰时间 / s | %.1f | %.1f | %.1f |\n',old.ice_onset_s(2),new.ice_onset_s(2),refs(2,3));
fprintf(fid,'| −25 ℃ | 峰值冰体积分数 | %.8f | %.8f | %.5f |\n\n',old.ice_volume_peak(2),new.ice_volume_peak(2),refs(2,4));
fprintf(fid,'## 数值验证与参数\n\n');
fprintf(fid,'优化退出标志为 %d；仅表示局部停止条件。\n\n',report.fit.exitflag);
fprintf(fid,'中→细网格检查通过：%d；细网格时间精度检查通过：%d。目标为最大电压差 5 mV、温度差 0.05 ℃、冰体积分数差 0.005。\n\n',report.meshPass,report.timePass);
fine=report.summary(report.summary.grid=="fine",:);
fprintf(fid,'在旧尺度 0.02 V / 0.5 ℃ 下，细网格加权误差为 %.6f；在本次尺度 0.02 V / 0.1 ℃ 下为 %.6f。不能跨权重直接比较优化目标值。\n\n',sum(fine.weighted_SSE),sum(fine.calibration_SSE));
fprintf(fid,'| 参数 | 数值 | 下界 | 上界 | 接近边界 |\n|---|---:|---:|---:|---:|\n');
p=report.parameters;
for j=1:height(p)
    fprintf(fid,'| %s | %.8g | %.8g | %.8g | %d |\n',p.names(j),p.values(j),p.lower(j),p.upper(j),p.at_bound(j));
end
fprintf(fid,'\n## 冰量假设与适用范围\n\n');
fprintf(fid,'起冰定义为每 0.2 s 输出的空间最大绝对孔隙冰体积分数首次超过 1e-6；峰值是时空最大 mi/rhoI，不是孔隙饱和度。截图的统计口径尚未核实。\n\n');
fprintf(fid,'基准成核阈值 0.03、冻结速率 1 s⁻¹ 均为假设。另在中网格比较阈值 0.02/0.03/0.04 与冻结速率 0.5/1/2 s⁻¹，其他参数固定；详情见 phase_assumption_sensitivity.csv。\n\n');
fprintf(fid,'电压与温度误差均是训练数据误差，尚无独立验证集。冰量没有实验观测，不能将更接近其他模型的结冰数值解释为实验验证通过。模型未保证每项 RMSE 都同时下降。\n');
end
