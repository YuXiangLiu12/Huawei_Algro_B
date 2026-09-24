function fit=q1_revise(root,mode,maxEvals)
% Revised local liquid/ice model. Reference-model ice values are never fitted.
if nargin<1,root=[];end
root=q1_project_root(root);
if nargin<2,mode='transport';end
if nargin<3,maxEvals=300;end
folder=fullfile(root,'results','q1_v06');if ~isfolder(folder),mkdir(folder);end
if strcmp(mode,'transport')
    previous=fullfile(root,'results','q1_v05','fit_transport_direct.mat');
    if isfile(previous)
        s=load(previous,'fit');c=s.fit.cfg;
    else
        c=q1_config();
        c.j0ref=.25;c.kineticHydrationExponent=2;
        c.boundWaterDiffusivityScale=2;c.fitTransport=true;
    end
    c.freezingClosure='delayed_nucleation';c.nucleationMode='local';
    c.liquidMigration=true;c.nucleationThreshold=.03;
    c.fitPhaseTransfer=false;active=1:5;
elseif strcmp(mode,'phase')
    s=load(fullfile(folder,'fit_transport.mat'),'fit');
    c=s.fit.cfg;c.fitPhaseTransfer=true;active=1:6;
else
    error('mode must be transport or phase');
end
% Working balance of the two observables, not measured uncertainty estimates.
c.fitResidualScales=[.02 .1];
c.runCalibration=true;c.version='q1-v0.6-local-phase';
c.calibrationCheckpoint=fullfile(folder,['checkpoint_' mode '.mat']);
if isfile(c.calibrationCheckpoint)
    p=load(c.calibrationCheckpoint,'checkpoint');
    c=q1_fit_config(p.checkpoint.theta,c);
elseif isfile(fullfile(folder,['fit_' mode '.mat']))
    p=load(fullfile(folder,['fit_' mode '.mat']),'fit');
    c=q1_fit_config(p.fit.theta,c);
end
% Always use the attachments at the current project location, not saved cases.
[cases,audit]=q1_read_data(root,c);
writetable(struct2table(audit),fullfile(folder,'input_audit.csv'));
fprintf('REVISED LOCAL PHASE FIT: %s; voltage scale .02 V, temperature scale .1 C\n',mode);
fit=q1_calibrate(cases,c,'coarse',maxEvals,1,active);
save(fullfile(folder,['fit_' mode '.mat']),'fit','cases');
disp(fit.detail);disp(fit.theta);disp(fit.resnorm);
end
