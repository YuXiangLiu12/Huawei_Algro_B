function q1_stage1(root)
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
folder=fullfile(root,'results','q1_v05');if ~isfolder(folder),mkdir(folder);end
c=q1_config();[cases,audit]=q1_read_data(root,c);
checks=q1_checks(root);writetable(checks,fullfile(folder,'module_checks.csv'));
kin=q1_fit_kinetics(cases,c);save(fullfile(folder,'kinetics.mat'),'kin');
c.j0ref=kin.j0_298;c.Ea=kin.Ea_Jmol;
fprintf('STAGE 1A: kinetics and conductivity, direct ice, beta=0\n');
fit=q1_calibrate(cases,c,'coarse',110,1,1:3);
save(fullfile(folder,'fit_stage1a.mat'),'fit','cases','audit');
disp(struct2table(fit.detail));
fprintf('STAGE 1B: allow existing hydration-dependent kinetics\n');
c=fit.cfg;c.kineticHydrationExponent=2;
fit=q1_calibrate(cases,c,'coarse',160,1,1:4);
save(fullfile(folder,'fit_stage1b.mat'),'fit','cases','audit');
disp(struct2table(fit.detail));
fprintf('STAGE 1 COMPLETE\n');
end
