function fit=q1_fit_transport(root,closure,maxEvals)
if nargin<1,root=[];end
root=q1_project_root(root);
if nargin<2,closure='direct';end
if nargin<3,maxEvals=600;end
folder=fullfile(root,'results','q1_v05');c=q1_config();[cases,~]=q1_read_data(root,c);
if ~isfolder(folder),mkdir(folder);end
c.j0ref=.25;c.Ea=67000;c.fCL=1;c.kineticHydrationExponent=2;
c.boundWaterDiffusivityScale=2;c.fitTransport=true;
if strcmp(closure,'local')
    c.freezingClosure='delayed_nucleation';c.nucleationMode='local';
    c.liquidMigration=true;c.nucleationThreshold=.03;
elseif ~strcmp(closure,'direct')
    error('closure must be direct or local');
end
fprintf('TRANSPORT JOINT FIT: %s\n',closure);
fitfile=fullfile(folder,['fit_transport_' closure '.mat']);
checkpointFile=fullfile(folder,['checkpoint_transport_' closure '.mat']);
if isfile(fitfile)
    previous=load(fitfile,'fit');c=previous.fit.cfg;
elseif isfile(checkpointFile)
    previous=load(checkpointFile,'checkpoint');c=previous.checkpoint.cfg;
end
c.calibrationCheckpoint=checkpointFile;
fit=q1_calibrate(cases,c,'coarse',maxEvals,1,1:5);
save(fitfile,'fit','cases');
disp(struct2table(fit.detail));disp(fit.theta);
end
