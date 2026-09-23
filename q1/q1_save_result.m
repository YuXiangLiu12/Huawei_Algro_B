function out=q1_save_result(cases,c,folder,label,level)
% Export only complete valid runs; declare fitted-vs-held-out status in label.
if nargin<5,level='medium';end
if ~isfolder(folder),mkdir(folder);end
c.runCalibration=true;c.version='q1-v0.5-calibrated';
out.cfg=c;out.cases=cases;
for k=1:numel(cases)
    run=q1_simulate(cases(k),c,level);out.runs(k)=run;
    assert(strcmp(run.status,'VALID'),'Invalid final trajectory: %s',run.event);
    d=cases(k);T=struct2table(run.obs);
    assert(height(T)==numel(d.t),'Incomplete final trajectory');
    T.V_exp=d.V_exp;T.T_exp_C=d.T_exp_C;
    T.j_Acm2=d.j_Acm2;
    T.water_residual=run.balance.waterResidual';
    T.energy_residual=run.balance.energyResidual';
    writetable(T,fullfile(folder,sprintf('case%d_observations.csv',k)));
    save(fullfile(folder,sprintf('case%d_run.mat',k)),'run','d','c','-v7.3');
end
status=table(string({out.runs.status})',string({out.runs.event})',[out.runs.t_end]', ...
    'VariableNames',{'status','event','t_end'});
writetable(status,fullfile(folder,'solver_status.csv'));
q1_export_results(out,folder,label);
q1_export_tables(out,folder,label);
rows=cell(0,13);
for k=1:numel(cases)
    d=cases(k);o=out.runs(k).obs;v=[o.V_sim]';tt=[o.T_seven_C]';
    ice=[o.ice_volume_max]';t=[o.time_s]';idx=find(ice>1e-6,1);
    onset=NaN;if ~isempty(idx),onset=t(idx);end
    early=d.t<=6;late=d.t>6;
    rows(end+1,:)={string(d.name),string(label), ...
        sqrt(mean((v-d.V_exp).^2)),sqrt(mean((tt-d.T_exp_C).^2)), ...
        sqrt(mean((v(early)-d.V_exp(early)).^2)), ...
        sqrt(mean((v(late)-d.V_exp(late)).^2)),onset,max(ice), ...
        max([o.ice_saturation_max]),max([o.liquid_saturation_max]), ...
        out.runs(k).balance.ew,out.runs(k).balance.eE,string(level)};
end
out.summary=cell2table(rows,'VariableNames',{'case_name','result_type','V_RMSE', ...
    'T_RMSE_K','V_RMSE_0_6','V_RMSE_6_end','ice_onset_s','ice_volume_peak', ...
    'ice_saturation_peak','liquid_saturation_peak','water_error','energy_error','grid'});
writetable(out.summary,fullfile(folder,'summary.csv'));
end
