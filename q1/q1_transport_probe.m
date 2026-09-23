function report=q1_transport_probe(root)
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
c=q1_config();[d,~]=q1_read_data(root,c);c.j0ref=.25;c.kineticHydrationExponent=2;
c.reltol=2e-5;c.abstol=2e-8;c.maxstep=.25;rows=cell(0,6);
for scale=[1 2 3]
    c.boundWaterDiffusivityScale=scale;
    for k=1:2
        r=q1_simulate(d(k),c,'coarse');
        rows(end+1,:)={scale,string(d(k).name),string(r.status), ...
            sqrt(mean(([r.obs.V_sim]'-d(k).V_exp).^2)), ...
            sqrt(mean(([r.obs.T_seven_C]'-d(k).T_exp_C).^2)),r.obs(end).lambda_aCL};
    end
end
report=cell2table(rows,'VariableNames',{'Db_scale','case_name','status','V_RMSE','T_RMSE_C','lambda_aCL_end'});
folder=fullfile(root,'results','q1_v05');if ~isfolder(folder),mkdir(folder);end
writetable(report,fullfile(folder,'transport_probe.csv'));disp(report);
end
