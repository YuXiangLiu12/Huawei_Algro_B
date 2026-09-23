function fit=q1_revise(root,mode)
% Revised local liquid/ice model. Reference-model ice values are never fitted.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
if nargin<2,mode='transport';end
folder=fullfile(root,'results','q1_v06');if ~isfolder(folder),mkdir(folder);end
if strcmp(mode,'transport')
    previous=fullfile(root,'results','q1_v05','fit_transport_direct.mat');
    if isfile(previous)
        s=load(previous,'fit','cases');c=s.fit.cfg;cases=s.cases;
    else
        c=q1_config();[cases,~]=q1_read_data(root,c);
        c.j0ref=.25;c.kineticHydrationExponent=2;
        c.boundWaterDiffusivityScale=2;c.fitTransport=true;
    end
    c.freezingClosure='delayed_nucleation';c.nucleationMode='local';
    c.liquidMigration=true;c.nucleationThreshold=.03;
    c.fitPhaseTransfer=false;active=1:5;
elseif strcmp(mode,'phase')
    s=load(fullfile(folder,'fit_transport.mat'),'fit','cases');
    c=s.fit.cfg;cases=s.cases;c.fitPhaseTransfer=true;active=1:6;
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
end
fprintf('REVISED LOCAL PHASE FIT: %s; voltage scale .02 V, temperature scale .1 C\n',mode);
fit=q1_calibrate(cases,c,'coarse',300,1,active);
save(fullfile(folder,['fit_' mode '.mat']),'fit','cases');
disp(fit.detail);disp(fit.theta);disp(fit.resnorm);
end
