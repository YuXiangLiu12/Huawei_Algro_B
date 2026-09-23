function cc=q1_fit_config(theta,c)
% Well-conditioned coordinates: log(j0 at -20 C), Ea/1e4, log(fCL), beta.
cc=c;
if numel(theta)==2 % Compatibility with pre-v0.5 diagnostic callers.
    cc.j0ref=exp(theta(1)); cc.fCL=exp(theta(2)); return;
end
cc.Ea=1e4*theta(2);
cc.j0ref=exp(theta(1))*exp(cc.Ea/cc.R*(1/253.15-1/298.15));
cc.fCL=exp(theta(3)); cc.kineticHydrationExponent=theta(4);
if numel(theta)>=5,cc.boundWaterDiffusivityScale=exp(theta(5));end
if numel(theta)>=6,cc.kbl=exp(theta(6));end
end
