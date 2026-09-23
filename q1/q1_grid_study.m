function summary=q1_grid_study(root,c,caseIndex)
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
if nargin<2,c=q1_config();end
if nargin<3,caseIndex=1:2;end
[cases,~]=q1_read_data(root,c);
rows=cell(0,10);
for k=caseIndex
    d=cases(k); prior=[]; priorName='';
    levels={'coarse','medium','fine'};
    for z=1:numel(levels)
        level=levels{z}; tic; run=q1_simulate(d,c,level); elapsed=toc;
        dv=NaN; dt=NaN; dice=NaN;
        if ~isempty(prior) && isfield(run,'obs') && ...
                numel(run.obs)==numel(prior.obs)
            dv=max(abs([run.obs.V_sim]-[prior.obs.V_sim]));
            dt=max(abs([run.obs.T_seven_C]-[prior.obs.T_seven_C]));
            dice=max(abs([run.obs.ice_volume_max]-[prior.obs.ice_volume_max]));
        end
        rows(end+1,:)={string(d.name),string(level),string(run.status), ...
            run.t_end,elapsed,run.balance.ew,run.balance.eE,dv,dt,dice}; %#ok<AGROW>
        prior=run; priorName=level; %#ok<NASGU>
    end
end
summary=cell2table(rows,'VariableNames',{'case_name','grid','status','t_end_s', ...
    'elapsed_s','water_error','energy_error','dV_from_coarser', ...
    'dT_from_coarser_K','dice_from_coarser'});
folder=fullfile(root,'results','checks'); if ~isfolder(folder),mkdir(folder);end
tag=sprintf('%d_',caseIndex); tag=tag(1:end-1);
writetable(summary,fullfile(folder,['grid_study_case' tag '.csv']));
end
