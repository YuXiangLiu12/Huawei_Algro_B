function q1_fit_transport(root,closure)
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
if nargin<2,closure='direct';end
folder=fullfile(root,'results','q1_v05');c=q1_config();[cases,~]=q1_read_data(root,c);
c.j0ref=.25;c.Ea=67000;c.fCL=1;c.kineticHydrationExponent=2;
c.boundWaterDiffusivityScale=2;c.fitTransport=true;
if strcmp(closure,'local')
    c.freezingClosure='delayed_nucleation';c.nucleationMode='local';
    c.liquidMigration=true;c.nucleationThreshold=.03;
elseif ~strcmp(closure,'direct')
    error('closure must be direct or local');
end
fprintf('TRANSPORT JOINT FIT: %s\n',closure);
fit=q1_calibrate(cases,c,'coarse',180,1,1:5);
save(fullfile(folder,['fit_transport_' closure '.mat']),'fit','cases');
disp(struct2table(fit.detail));disp(fit.theta);
end
