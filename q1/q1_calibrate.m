function fit=q1_calibrate(cases,c,level,maxEvals,nStarts,active)
% Joint voltage-temperature calibration. Both cases share one parameter set.
if nargin<3,level='coarse';end
if nargin<4,maxEvals=160;end
if nargin<5,nStarts=1;end
if nargin<6,active=1:4;end
lb=[log(1e-5);2;log(.1);0]; ub=[log(.1);9;log(5);6];
xbase=[log(c.j0ref)-c.Ea/c.R*(1/253.15-1/298.15); ...
    c.Ea/1e4;log(c.fCL);c.kineticHydrationExponent];
if isfield(c,'fitTransport') && c.fitTransport
    lb=[lb;log(.3)];ub=[ub;log(5)];
    xbase=[xbase;log(c.boundWaterDiffusivityScale)];
end
if isfield(c,'fitPhaseTransfer') && c.fitPhaseTransfer
    assert(numel(xbase)==5,'Phase transfer calibration requires fitTransport');
    lb=[lb;log(.02)];ub=[ub;log(5)];xbase=[xbase;log(c.kbl)];
end
xbase=min(ub,max(lb,xbase));
% Keep expensive calibration tolerances separate from the verification run.
cc=c; cc.reltol=2e-5; cc.abstol=2e-8; cc.maxstep=.25;
fit.attempts=repmat(struct(),1,nStarts);
history=struct('start',{},'iteration',{},'funccount',{},'resnorm',{},'theta',{});
for i=1:nStarts
    x0=xbase;
    if i>1,x0(3)=min(ub(3),max(lb(3),x0(3)+(-1)^i*.6));end
    pack=@(z) assemble(z,xbase,active);
    objective=@(z) q1_residual(pack(z),cases,cc,level);
    if license('test','optimization_toolbox') && exist('lsqnonlin','file')==2
        opt=optimoptions('lsqnonlin','Display','iter','MaxFunctionEvaluations',maxEvals, ...
            'FiniteDifferenceStepSize',.002,'FunctionTolerance',2e-5, ...
            'StepTolerance',2e-4,'OptimalityTolerance',1e-3);
        opt.OutputFcn=@saveProgress;
        [z,resnorm,residual,flag,out,~,jac]=lsqnonlin(objective,x0(active),lb(active),ub(active),opt);
    else
        opt=optimset('Display','iter','MaxFunEvals',maxEvals,'TolX',2e-4,'TolFun',2e-5);
        trans=@(z) lb(active)+(ub(active)-lb(active))./(1+exp(-z));
        opt.OutputFcn=@(zz,values,state)saveProgress(trans(zz),values,state);
        ratio=(x0(active)-lb(active))./(ub(active)-lb(active));
        ratio=min(1-1e-5,max(1e-5,ratio));
        z0=log(ratio./(1-ratio));
        [zz,resnorm,flag,out]=fminsearch(@(zz)sum(objective(trans(zz)).^2),z0,opt);
        z=trans(zz); residual=objective(z); jac=[];
    end
    x=pack(z);
    fit.attempts(i).x=x; fit.attempts(i).resnorm=resnorm;
    fit.attempts(i).exitflag=flag; fit.attempts(i).output=out;
    fit.attempts(i).residual=residual; fit.attempts(i).jacobian=jac;
end
[~,best]=min([fit.attempts.resnorm]); fit.theta=fit.attempts(best).x;
fit.cfg=q1_fit_config(fit.theta,c); fit.j0ref=fit.cfg.j0ref;
fit.j0_253=exp(fit.theta(1));fit.fCL=fit.cfg.fCL;
fit.Ea_Jmol=fit.cfg.Ea;fit.beta=fit.cfg.kineticHydrationExponent;
fit.Db_scale=1;
if isfield(fit.cfg,'boundWaterDiffusivityScale'),fit.Db_scale=fit.cfg.boundWaterDiffusivityScale;end
fit.resnorm=fit.attempts(best).resnorm; fit.exitflag=fit.attempts(best).exitflag;
fit.level=level;fit.active=active;fit.lower=lb;fit.upper=ub;
fit.trainingCases=string({cases.name});fit.weights=[.02 .5];
if isfield(c,'fitResidualScales'),fit.weights=c.fitResidualScales;end
fit.kbl=fit.cfg.kbl;
fit.note='Joint calibration errors, not independent experimental validation; ice not fitted.';
fit.history=history;fit.converged=fit.exitflag>0;
[~,fit.detail]=q1_residual(fit.theta,cases,cc,level);
    function stop=saveProgress(z,values,state)
        stop=false;
        if ~any(strcmp(state,{'init','iter','done'})),return;end
        checkpoint.theta=assemble(z,xbase,active);
        checkpoint.cfg=q1_fit_config(checkpoint.theta,c);
        checkpoint.start=i;checkpoint.state=state;checkpoint.level=level;
        checkpoint.optimValues=values;
        value=NaN;
        if isfield(values,'resnorm'),value=values.resnorm;
        elseif isfield(values,'fval'),value=values.fval;end
        history(end+1)=struct('start',i,'iteration',values.iteration, ...
            'funccount',values.funccount,'resnorm',value,'theta',checkpoint.theta);
        checkpoint.history=history;
        if isfield(c,'calibrationCheckpoint') && ~isempty(c.calibrationCheckpoint)
            save(c.calibrationCheckpoint,'checkpoint');
        end
    end
end
function x=assemble(z,base,active)
x=base;x(active)=z;
end
