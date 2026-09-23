function [r,detail]=q1_residual(theta,cases,c,level)
% See q1_fit_config for fit coordinates. Failed runs get a fixed-length penalty.
if nargin<4,level='coarse';end
cc=q1_fit_config(theta,c);
scales=[.02 .5];
if isfield(c,'fitResidualScales'),scales=c.fitResidualScales;end
assert(numel(scales)==2 && all(isfinite(scales)) && all(scales>0), ...
    'fitResidualScales must contain positive voltage and temperature scales');
r=[]; detail=repmat(struct('status','','event','','V_RMSE',NaN,'T_RMSE_C',NaN),1,numel(cases));
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
        rv=(V-d.V_exp)/scales(1); rt=(T-d.T_exp_C)/scales(2);
        detail(k).V_RMSE=sqrt(mean((V-d.V_exp).^2));
        detail(k).T_RMSE_C=sqrt(mean((T-d.T_exp_C).^2));
    else
        remaining=max(0,1-run.t_end/d.t(end));
        if isnan(remaining),remaining=1;end
        rv=50*(1+remaining)*ones(n,1); rt=rv;
    end
    r=[r;rv;rt]; %#ok<AGROW>
end
end
