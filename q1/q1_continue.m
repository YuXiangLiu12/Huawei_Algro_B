function out=q1_continue(root,maxEvals)
% Resume saved local fits, verify both parameter counts, then compare equally.
% maxEvals=0 verifies saved fits without running either optimizer.
if nargin<1,root=[];end
root=q1_project_root(root);
if nargin<2,maxEvals=300;end
base=fullfile(root,'results','q1_v06');
if ~isfolder(base),mkdir(base);end
modes={'transport','phase'};
for z=1:2
    mode=modes{z};file=fullfile(base,['fit_' mode '.mat']);
    if maxEvals>0
        done=false;
        if isfile(file)
            s=load(file,'fit','cases');
            [current,~]=q1_read_data(root,s.fit.cfg);
            same=true;
            for k=1:numel(current)
                fields={'t','I_A','V_exp','T_exp_C','j_SI'};
                for j=1:numel(fields)
                    same=same && isequal(current(k).(fields{j}),s.cases(k).(fields{j}));
                end
            end
            done=s.fit.exitflag>0 && same;
        end
        if ~done,q1_revise(root,mode,maxEvals);end
    end
    assert(isfile(file),'Missing %s. Run q1_main([],''revised'') first.',file);
    out.(mode)=q1_verify_transport(root,['local_' mode],file, ...
        fullfile(base,['verification_' mode]));
end
% Rerun the saved direct model on the SAME fine mesh and tolerances.
previous=load(fullfile(root,'results','q1_v05','fit_transport_direct.mat'),'fit');
[cases,~]=q1_read_data(root,previous.fit.cfg);
baseline=q1_save_result(cases,previous.fit.cfg,fullfile(base,'baseline_direct_fine'), ...
    'previous_direct_fixed_parameters','fine');
out.baseline=baseline.summary;
out.comparison=q1_compare_revision(root,out);
% The assumptions are not calibration targets. Check the selected candidate.
selected=out.comparison.selectedMode;
out.(selected)=q1_revision_report(root,selected,out.(selected));
out.comparison=q1_compare_revision(root,out); % Include completed sensitivity results.
save(fullfile(base,'comparison.mat'),'out');
end
