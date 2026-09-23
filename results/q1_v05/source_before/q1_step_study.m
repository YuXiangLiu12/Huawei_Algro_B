function summary=q1_step_study(root,caseIndex)
if nargin<1,root=fileparts(mfilename('fullpath'));end
if nargin<2,caseIndex=1;end
c=q1_config(); [cases,~]=q1_read_data(root,c); d=cases(caseIndex);
steps=[.05 .025 .0125 .0125]; tolerances=[1e-6 1e-6 1e-6 1e-7];
rows=cell(numel(steps),9); reference=[];
for k=1:numel(steps)
    cc=c; cc.maxstep=steps(k); cc.reltol=tolerances(k);
    run=q1_simulate(d,cc,'medium');
    dv=NaN; dt=NaN; dice=NaN;
    if ~isempty(reference) && numel(run.obs)==numel(reference.obs)
        dv=max(abs([run.obs.V_sim]-[reference.obs.V_sim]));
        dt=max(abs([run.obs.T_seven_C]-[reference.obs.T_seven_C]));
        dice=max(abs([run.obs.ice_volume_max]-[reference.obs.ice_volume_max]));
    end
    rows(k,:)={string(d.name),steps(k),tolerances(k),string(run.status), ...
        run.balance.ew,run.balance.eE,dv,dt,dice};
    reference=run;
end
summary=cell2table(rows,'VariableNames',{'case_name','max_step_s','reltol', ...
    'status','water_error','energy_error','dV_previous','dT_previous_K', ...
    'dice_previous'});
folder=fullfile(root,'results','checks'); if ~isfolder(folder),mkdir(folder);end
writetable(summary,fullfile(folder,sprintf('step_study_case%d.csv',caseIndex)));
end
