function [r,detail]=q1_residual(theta,cases,c,level)
% See q1_fit_config for fit coordinates. Failed runs get a fixed-length penalty.
if nargin<4,level='coarse';end
cc=q1_fit_config(theta,c);
r=[]; detail=repmat(struct('status','','event','','V_RMSE',NaN,'T_RMSE',NaN),1,numel(cases));
for k=1:numel(cases)
    d=cases(k); run=q1_simulate(d,cc,level); n=numel(d.t);
    detail(k).status=run.status; detail(k).event=run.event;
    if strcmp(run.status,'VALID') && isfield(run,'obs') && numel(run.obs)==n
        V=[run.obs.V_sim]';
        if strcmp(cc.temperatureAverage,'seven_layer')
            T=[run.obs.T_seven_C]';
        else
            T=[run.obs.T_MEA_C]';
        end
        rv=(V-d.V_exp)/.02; rt=(T-d.T_exp_C)/.5;
        detail(k).V_RMSE=sqrt(mean((V-d.V_exp).^2));
        detail(k).T_RMSE=sqrt(mean((T-d.T_exp_C).^2));
    else
        remaining=max(0,1-run.t_end/d.t(end));
        if isnan(remaining),remaining=1;end
        rv=50*(1+remaining)*ones(n,1); rt=rv;
    end
    r=[r;rv;rt]; %#ok<AGROW>
end
end
