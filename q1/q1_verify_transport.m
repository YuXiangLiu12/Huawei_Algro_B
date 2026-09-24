function report=q1_verify_transport(root,closure,fitfile,folder)
% Verify one saved fit without requiring historical fits or source snapshots.
% Thresholds below are numerical acceptance targets, not experimental errors.
if nargin<1,root=[];end
root=q1_project_root(root);
if nargin<2,closure='direct';end
base=fullfile(root,'results','q1_v05');
if nargin<3,fitfile=fullfile(base,['fit_transport_' closure '.mat']);end
s=load(fitfile,'fit','cases');
fit=s.fit;c=fit.cfg;[cases,~]=q1_read_data(root,c);
if isfield(c,'calibrationCheckpoint'),c.calibrationCheckpoint='';end
if nargin<4,folder=fullfile(base,['verification_' closure]);end
if ~isfolder(folder),mkdir(folder);end
limits=[.005 .05 .005]; % V, K, absolute ice volume fraction
levels={'coarse','medium','fine'};rows=cell(0,15);comparisons=cell(0,8);
for z=1:numel(levels)
    level=levels{z};
    fprintf('VERIFY %s grid at fixed fitted parameters\n',level);
    out=q1_save_result(cases,c,fullfile(folder,level), ...
        ['transport_' closure '_joint_calibration'],level);
    runs{z}=out.runs; %#ok<AGROW>
    if z==3,final=out;end
    for k=1:numel(cases)
        run=out.runs(k);d=cases(k);o=run.obs;
        rv=([o.V_sim]'-d.V_exp)/.02;
        rt=([o.T_seven_C]'-d.T_exp_C)/.5;
        rows(end+1,:)={string(d.name),string(level),run.grid.N,run.index.n, ...
            string(run.status),sum(rv.^2)+sum(rt.^2), ...
            sqrt(mean((rv*.02).^2)),sqrt(mean((rt*.5).^2)), ...
            run.balance.ew,run.balance.eE,min([o.lambda_aCL]), ...
            o(end).lambda_aCL,max([o.ice_volume_max]),run.t_end, ...
            sum((rv*.02/fit.weights(1)).^2)+sum((rt*.5/fit.weights(2)).^2)}; %#ok<AGROW>
        if z>1
            comparisons(end+1,:)=compare(run,runs{z-1}(k),d, ...
                [level '_vs_' levels{z-1}],limits); %#ok<AGROW>
        end
    end
    summary=cell2table(rows,'VariableNames',{'case_name','grid','cells','states', ...
        'status','weighted_SSE','V_RMSE','T_RMSE_C','water_error','energy_error', ...
        'lambda_aCL_min','lambda_aCL_end','ice_volume_peak','t_end_s','calibration_SSE'});
    writetable(summary,fullfile(folder,'grid_summary.csv'));
end
% Check time accuracy on the fine mesh used in the delivered trajectories.
cc=c;cc.maxstep=.025;cc.reltol=2e-7;cc.abstol=2e-10;
for k=1:numel(cases)
    fprintf('VERIFY tighter time integration: case %d\n',k);
    run=q1_simulate(cases(k),cc,'fine');
    assert(strcmp(run.status,'VALID'),'Tighter integration failed: %s',run.event);
    comparisons(end+1,:)=compare(run,runs{3}(k),cases(k),'tight_vs_fine',limits);
    save(fullfile(folder,sprintf('case%d_tight.mat',k)),'run','cc','-v7.3');
end
checks=cell2table(comparisons,'VariableNames',{'case_name','comparison', ...
    'max_dV','max_dT_K','max_dice','water_error','energy_error','pass'});
writetable(checks,fullfile(folder,'numerical_checks.csv'));
names=["j0_253_A_m2";"Ea_J_mol";"fCL";"beta";"Db_scale"];
values=[fit.j0_253;fit.Ea_Jmol;fit.fCL;fit.beta;fit.Db_scale];
lower=[exp(fit.lower(1));fit.lower(2)*1e4;exp(fit.lower(3)); ...
    fit.lower(4);exp(fit.lower(5))];
upper=[exp(fit.upper(1));fit.upper(2)*1e4;exp(fit.upper(3)); ...
    fit.upper(4);exp(fit.upper(5))];
if numel(fit.theta)>=6
    names=[names;"bound_to_liquid_rate_s_inv"];values=[values;fit.cfg.kbl];
    lower=[lower;exp(fit.lower(6))];upper=[upper;exp(fit.upper(6))];
end
at_bound=min(fit.theta-fit.lower,fit.upper-fit.theta)<1e-3;
parameters=table(names,values,lower,upper,at_bound);
writetable(parameters,fullfile(folder,'fitted_parameters.csv'));
if isfield(fit,'history') && ~isempty(fit.history)
    h=fit.history;history=table([h.start]',[h.iteration]',[h.funccount]',[h.resnorm]', ...
        'VariableNames',{'start','iteration','funccount','weighted_SSE'});
    writetable(history,fullfile(folder,'optimization_history.csv'));
end
report.fit=fit;report.summary=summary;report.checks=checks;report.parameters=parameters;
report.numericalLimits=limits;
report.optimizerStoppedSuccessfully=fit.exitflag>0;
report.meshPass=all(checks.pass(checks.comparison=="fine_vs_medium"));
report.timePass=all(checks.pass(checks.comparison=="tight_vs_fine"));
save(fullfile(folder,'verification.mat'),'report');
q1_plot_mesh(final.cases,runs,folder);
fid=fopen(fullfile(folder,'verification_report.md'),'w','n','UTF-8');
cleanup=onCleanup(@()fclose(fid));
fprintf(fid,'# Direct verification of the %s ice model\n\n',closure);
fprintf(fid,'Joint fit to both temperatures; these are training errors, not held-out validation. Ice is not observed.\n\n');
fprintf(fid,'Optimizer exit flag: %d. Coarse calibration weighted SSE: %.9g.\n\n',fit.exitflag,fit.resnorm);
fprintf(fid,'Calibration residual scales: %.6g V and %.6g deg C.\n\n',fit.weights);
[~,best]=min([fit.attempts.resnorm]);
fprintf(fid,'Optimizer message: %s\n\n',fit.attempts(best).output.message);
fprintf(fid,'Fixed-parameter, strict-tolerance SSE with common scales 0.02 V and 0.5 deg C (for cross-model comparison):\n\n');
for z=1:3
    fprintf(fid,'- %s: %.9g\n',levels{z},sum(summary.weighted_SSE(summary.grid==string(levels{z}))));
end
fprintf(fid,'\nNumerical targets (maximum over observed times): %.3g V, %.3g K, %.3g absolute ice volume fraction.\n\n',limits);
fprintf(fid,'Fine-versus-medium mesh check passed: %d. Fine-grid time check passed: %d.\n\n',report.meshPass,report.timePass);
fprintf(fid,'Parameters near a calibration bound:\n\n');
for j=find(at_bound)'
    fprintf(fid,'- %s = %.9g, allowed range [%.9g, %.9g]\n', ...
        names(j),values(j),lower(j),upper(j));
end
fprintf(fid,'\nBound-active parameters should not be interpreted as independently identified physical constants.\n\n');
fprintf(fid,'A positive optimizer exit flag is a local stopping condition, not evidence of a global optimum. If the mesh check fails, the fit is not mesh independent and needs finer-grid calibration.\n');
disp(parameters);disp(summary);disp(checks);
end

function row=compare(a,b,d,label,limits)
assert(numel(a.obs)==numel(d.t) && numel(b.obs)==numel(d.t),'Incomplete trajectory');
delta=[max(abs([a.obs.V_sim]-[b.obs.V_sim])), ...
    max(abs([a.obs.T_seven_C]-[b.obs.T_seven_C])), ...
    max(abs([a.obs.ice_volume_max]-[b.obs.ice_volume_max]))];
row={string(d.name),string(label),delta(1),delta(2),delta(3), ...
    a.balance.ew,a.balance.eE,all(delta<=limits)};
end
