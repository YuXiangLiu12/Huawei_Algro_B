function fit=q1_calibrate(cases,c,level,maxEvals,nStarts)
if nargin<3,level='coarse';end
if nargin<4,maxEvals=40;end
if nargin<5,nStarts=3;end
lb=log([1e-4;.3]); ub=log([10;3]);
starts=[log([.25;1]),log([.02;.6]),log([1;2])];
starts=starts(:,1:min(nStarts,size(starts,2)));
fit.attempts=repmat(struct(),1,size(starts,2));
for i=1:size(starts,2)
    x0=min(ub,max(lb,starts(:,i)));
    if license('test','optimization_toolbox') && exist('lsqnonlin','file')==2
        opt=optimoptions('lsqnonlin','Display','off','MaxFunctionEvaluations',maxEvals, ...
            'FunctionTolerance',1e-5,'StepTolerance',1e-5);
        [x,resnorm,~,flag,out]=lsqnonlin(@(z)q1_residual(z,cases,c,level),x0,lb,ub,opt);
    else
        opt=optimset('Display','off','MaxFunEvals',maxEvals,'TolX',1e-4);
        trans=@(z) lb+(ub-lb)./(1+exp(-z));
        z0=log((x0-lb+1e-6)./(ub-x0+1e-6));
        [z,resnorm,flag,out]=fminsearch(@(z)sum(q1_residual(trans(z),cases,c,level).^2),z0,opt);
        x=trans(z);
    end
    fit.attempts(i).x=x; fit.attempts(i).resnorm=resnorm;
    fit.attempts(i).exitflag=flag; fit.attempts(i).output=out;
end
[~,best]=min([fit.attempts.resnorm]);
fit.theta=fit.attempts(best).x; fit.j0ref=exp(fit.theta(1));
fit.fCL=exp(fit.theta(2)); fit.resnorm=fit.attempts(best).resnorm;
fit.level=level; fit.Ea_Jmol=c.Ea; fit.trainingCases=string({cases.name});
[~,fit.detail]=q1_residual(fit.theta,cases,c,level);
end
