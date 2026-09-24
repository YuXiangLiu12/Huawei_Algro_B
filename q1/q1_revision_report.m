function report=q1_revision_report(root,mode,report)
% Verify a local-phase fit and test unobserved phase assumptions separately.
if nargin<1,root=[];end
root=q1_project_root(root);
if nargin<2,mode='phase';end
base=fullfile(root,'results','q1_v06');
fitfile=fullfile(base,['fit_' mode '.mat']);
folder=fullfile(base,['verification_' mode]);
if nargin<3
    report=q1_verify_transport(root,['local_' mode],fitfile,folder);
end
[cases,~]=q1_read_data(root,report.fit.cfg);rows={};
settings=[.02 1;.03 .5;.03 1;.03 2;.04 1];
for z=1:size(settings,1)
    c=report.fit.cfg;c.nucleationThreshold=settings(z,1);c.kfreeze=settings(z,2);
    for k=1:numel(cases)
        fprintf('SENSITIVITY %s: threshold %.3g, freezing %.3g, case %d\n', ...
            mode,c.nucleationThreshold,c.kfreeze,k);
        r=q1_simulate(cases(k),c,'medium');
        assert(strcmp(r.status,'VALID'),'Sensitivity run failed: %s',r.event);
        o=r.obs;ice=[o.ice_volume_max];idx=find(ice>1e-6,1);onset=NaN;
        if ~isempty(idx),onset=o(idx).time_s;end
        v=[o.V_sim]'-cases(k).V_exp;t=[o.T_seven_C]'-cases(k).T_exp_C;
        rows(end+1,:)={c.nucleationThreshold,c.kfreeze,string(cases(k).name), ...
            sqrt(mean(v.^2)),sqrt(mean(t.^2)),onset,max(ice), ...
            sum((v/.02).^2+(t/.5).^2),sum((v/.02).^2+(t/.1).^2)}; %#ok<AGROW>
    end
end
report.sensitivity=cell2table(rows,'VariableNames',{'nucleation_threshold', ...
    'kfreeze_s_inv','case_name','V_RMSE','T_RMSE_C','ice_onset_s','ice_peak', ...
    'weighted_SSE','calibration_SSE'});
writetable(report.sensitivity,fullfile(folder,'phase_assumption_sensitivity.csv'));
save(fullfile(folder,'revision_report.mat'),'report');
end
