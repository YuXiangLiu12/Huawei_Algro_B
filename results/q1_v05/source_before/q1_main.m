function out=q1_main(root)
if nargin<1,root=fileparts(mfilename('fullpath'));end
c=q1_config(); [cases,audit]=q1_read_data(root,c);
folder=fullfile(root,'results','q1_v04');
if ~isfolder(folder),mkdir(folder);end
auditTable=struct2table(audit); writetable(auditTable,fullfile(folder,'input_audit.csv'));
fit=q1_fit_kinetics(cases,c);
save(fullfile(folder,'kinetics.mat'),'fit');
out.cfg=c; out.cases=cases; out.audit=audit; out.kinetics=fit;
for k=1:2
    run=q1_simulate(cases(k),c,'medium'); out.runs(k)=run;
    if isfield(run,'obs')
        d=cases(k); T=struct2table(run.obs); n=height(T);
        T.V_exp=d.V_exp(1:n); T.T_exp_C=d.T_exp_C(1:n);
        T.water_residual=run.balance.waterResidual';
        T.energy_residual=run.balance.energyResidual';
        writetable(T,fullfile(folder,sprintf('case%d_observations.csv',k)));
        save(fullfile(folder,sprintf('case%d_run.mat',k)),'run','d','c','-v7.3');
    end
end
status=table(string({out.runs.status})',string({out.runs.event})', ...
    [out.runs.t_end]','VariableNames',{'status','event','t_end'});
writetable(status,fullfile(folder,'solver_status.csv'));
q1_export_results(out,folder,'unfitted_preview');
end
