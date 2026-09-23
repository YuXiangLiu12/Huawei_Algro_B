function p=q1_properties(s,c,g)
% Derived quantities are evaluated on full thermal grid; inactive cells remain zero.
N=g.N; p.eg=zeros(N,1); p.sl=zeros(N,1); p.si=zeros(N,1);
p.eg(g.porous)=g.eps(g.porous)-s.ml(g.porous)/c.rhoL-s.mi(g.porous)/c.rhoI;
p.sl(g.porous)=s.ml(g.porous)./(c.rhoL*g.eps(g.porous));
p.si(g.porous)=s.mi(g.porous)./(c.rhoI*g.eps(g.porous));
egEval=max(p.eg,1e-12); T=max(s.T,200);
p.cH=s.qH./egEval; p.cO=s.qO./egEval;
p.pv=s.mv*c.R.*T./(c.Mw*egEval);
p.rv=s.mv./egEval;
p.lambda=zeros(N,1); p.lambda(g.ion)=s.bn(g.ion)./(g.omega(g.ion)*c.Cm);
Tc=T-273.15;
p.psl=611.21*exp((18.678-Tc/234.5).*Tc./(257.14+Tc));
p.psi=611.15*exp((23.036-Tc/333.7).*Tc./(279.82+Tc));
ps=p.psl;
if c.liquidRoute || strcmp(c.freezingClosure,'delayed_nucleation') || ...
        strcmp(c.activitySurface,'liquid_all')
    % Liquid saturation pressure is used for the activity branch.
elseif strcmp(c.activitySurface,'ice_below_zero')
    ps(T<273.15)=p.psi(T<273.15);
else
    error('Unknown activity surface');
end
a=p.pv./ps+2*p.sl; p.activity=a; ae=min(3,max(0,a));
p.lambdaEq=14+1.4*(ae-1); m=ae<=1;
p.lambdaEq(m)=.043+17.81*ae(m)-39.85*ae(m).^2+36*ae(m).^3;
p.lambdaSat=inf(N,1); m=T<223.15; p.lambdaSat(m)=4.837;
m=T>=223.15 & T<273.15;
p.lambdaSat(m)=1./(-1.304+.01479*T(m)-3.594e-5*T(m).^2);
lam=p.lambda; poly=2.563-.33*lam+.0264*lam.^2-.000671*lam.^3;
p.Db=1e-10*exp(2416*(1/303.15-1./T)).*poly;
p.kappa=(.5139*lam-.326).*exp(1268*(1/303.15-1./T));
if strcmp(c.conductivity,'cold2222')
    p.kappa=100*exp(2222*(1/303-1./T)).*(.005139*lam-.00326);
end
p.kappa(g.ion)=p.kappa(g.ion).* (1+(c.fCL*g.omega(g.ion).^1.5-1).*(g.layer(g.ion)~=4));
p.kappa(~g.ion)=0;
% Reference binary gas diffusivities are explicit working assumptions (m^2/s).
Dscale=(T/298.15).^1.75.*max(p.eg,0).^1.5;
if ~c.iceFeedback
    Dscale=(T/298.15).^1.75.*g.eps.^1.5;
end
p.DH=c.DHref*Dscale; p.DO=c.DOref*Dscale; p.Dv=c.Dvref*Dscale;
% Volume-fraction weighted k; the dry layer value remains the reference.
p.k=g.k; m=g.porous;
gasFraction=max(p.eg(m),0); liquidFraction=max(s.ml(m)/c.rhoL,0);
iceFraction=max(s.mi(m)/c.rhoI,0);
ksolid=(g.k(m)-g.eps(m)*c.kGas)./(1-g.eps(m));
p.k(m)=(1-g.eps(m)).*ksolid+gasFraction*c.kGas+ ...
    liquidFraction*c.kLiquid+iceFraction*c.kIce;
end
