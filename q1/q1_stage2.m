function q1_stage2(root)
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
folder=fullfile(root,'results','q1_v05');
s=load(fullfile(folder,'fit_stage1a.mat'));cases=s.cases;
k=load(fullfile(folder,'kinetics.mat'));c=q1_config();
c.Ea=k.kin.Ea_Jmol;c.j0ref=k.kin.j0_298;c.kineticHydrationExponent=2;
c.freezingClosure='delayed_nucleation';c.nucleationMode='local';
c.liquidMigration=true;c.nucleationThreshold=.03;
fprintf('STAGE 2: supercooled water, local nucleation, transport; threshold is assumed\n');
fit=q1_calibrate(cases,c,'coarse',150,1,1:4);
save(fullfile(folder,'fit_stage2.mat'),'fit','cases');
disp(struct2table(fit.detail));
end
