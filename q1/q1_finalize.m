function out=q1_finalize(root)
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
folder=fullfile(root,'results','q1_v05');
a=load(fullfile(folder,'fit_stage1a.mat'));b=load(fullfile(folder,'fit_stage1b.mat'));
z=load(fullfile(folder,'fit_stage2.mat'));cases=b.cases;
u=load(fullfile(folder,'fit_transport_direct.mat'));
v=load(fullfile(folder,'fit_transport_local.mat'));
fits={a.fit,b.fit,z.fit,u.fit,v.fit};
labels={'stage1a_baseline','stage1b_direct','stage2_local','transport_direct','transport_local'};
% All comparisons use the same data and residual scaling. Complexity penalty
% is 2 per fitted parameter, so an extra parameter must improve the objective.
scores=cellfun(@(f)f.resnorm+2*numel(f.active),fits);
[~,best]=min(scores);selected=fits{best};selectedStage=labels{best};
selected.cfg.runCalibration=true;selected.cfg.version='q1-v0.5-calibrated';
allSummary=table();
for i=1:numel(fits)
    fprintf('VERIFY %s on medium grid\n',labels{i});
    result=q1_save_result(cases,fits{i}.cfg,fullfile(folder,labels{i}), ...
        [labels{i} '_joint_calibration'],'medium');
    allSummary=[allSummary;result.summary];
    if strcmp(labels{i},selectedStage),out=result;end
end
writetable(allSummary,fullfile(folder,'stage_comparison.csv'));
% Final result remains a declared joint fit, not independent validation.
finaldir=fullfile(folder,'final');if ~isfolder(finaldir),mkdir(finaldir);end
srcdir=fullfile(folder,selectedStage);
files=dir(srcdir);
for i=1:numel(files)
    if ~files(i).isdir,copyfile(fullfile(srcdir,files(i).name),finaldir);end
end
q1_export_comparison(root,out,finaldir);
fit=selected;save(fullfile(finaldir,'selected_fit.mat'),'fit','selectedStage');
% Two-temperature mesh and time integration checks at the fitted parameters.
rows=cell(0,9);
for k=1:2
    base=out.runs(k);
    fine=q1_simulate(cases(k),selected.cfg,'fine');
    assert(strcmp(fine.status,'VALID'),'Fine grid run failed');
    rows(end+1,:)={string(cases(k).name),"fine_vs_medium", ...
        max(abs([fine.obs.V_sim]-[base.obs.V_sim])), ...
        max(abs([fine.obs.T_seven_C]-[base.obs.T_seven_C])), ...
        max(abs([fine.obs.ice_volume_max]-[base.obs.ice_volume_max])), ...
        sqrt(mean(([fine.obs.V_sim]'-cases(k).V_exp).^2)), ...
        sqrt(mean(([fine.obs.T_seven_C]'-cases(k).T_exp_C).^2)), ...
        fine.balance.ew,fine.balance.eE};
    cc=selected.cfg;cc.maxstep=.025;cc.reltol=2e-7;cc.abstol=2e-10;
    tight=q1_simulate(cases(k),cc,'medium');
    assert(strcmp(tight.status,'VALID'),'Time accuracy run failed');
    rows(end+1,:)={string(cases(k).name),"tighter_vs_medium", ...
        max(abs([tight.obs.V_sim]-[base.obs.V_sim])), ...
        max(abs([tight.obs.T_seven_C]-[base.obs.T_seven_C])), ...
        max(abs([tight.obs.ice_volume_max]-[base.obs.ice_volume_max])), ...
        sqrt(mean(([tight.obs.V_sim]'-cases(k).V_exp).^2)), ...
        sqrt(mean(([tight.obs.T_seven_C]'-cases(k).T_exp_C).^2)), ...
        tight.balance.ew,tight.balance.eE};
end
checks=cell2table(rows,'VariableNames',{'case_name','comparison','max_dV', ...
    'max_dT_K','max_dice','V_RMSE','T_RMSE_C','water_error','energy_error'});
writetable(checks,fullfile(folder,'numerical_checks.csv'));
q1_revision_checks(root);q1_local_phase_checks(root);
% Ice observations do not exist: vary the assumed threshold rather than fit it.
rows=cell(0,7);
for threshold=[.02 .03 .04]
    cc=v.fit.cfg;cc.nucleationThreshold=threshold;
    for k=1:2
        run=q1_simulate(cases(k),cc,'coarse');assert(strcmp(run.status,'VALID'));
        o=run.obs;ice=[o.ice_volume_max];index=find(ice>1e-6,1);onset=NaN;
        if ~isempty(index),onset=o(index).time_s;end
        rows(end+1,:)={threshold,string(cases(k).name),onset,max(ice), ...
            sqrt(mean(([o.V_sim]'-cases(k).V_exp).^2)), ...
            sqrt(mean(([o.T_seven_C]'-cases(k).T_exp_C).^2)),max([o.liquid_saturation_max])};
    end
end
sensitivity=cell2table(rows,'VariableNames',{'assumed_threshold','case_name','ice_onset_s', ...
    'ice_volume_peak','V_RMSE','T_RMSE_C','liquid_saturation_peak'});
writetable(sensitivity,fullfile(folder,'ice_assumption_sensitivity.csv'));
disp(allSummary);disp(checks);disp(sensitivity);
fprintf('SELECTED %s\n',selectedStage);
end
