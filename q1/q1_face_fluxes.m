function f=q1_face_fluxes(s,p,j,d,c,g,v)
% Every internal interface is evaluated once; rightward flux is positive.
N=g.N; f.T=zeros(N+1,1); f.H=f.T; f.O=f.T; f.v=f.T; f.b=f.T; f.l=f.T;
if nargin<7,gammaFace=g.gammaFace;else,gammaFace=v.gammaFace;end
if c.liquidMigration
    sl=min(1,max(p.sl,0));
    K0=c.K0Layer(g.layer)'; angle=c.thetaLayer(g.layer)';
    pc=c.surfaceTension*cosd(angle).*sqrt(g.eps./max(K0,realmin)).* ...
        (1.417*sl-2.120*sl.^2+1.263*sl.^3);
    mob=c.rhoL*K0.*sl.^3/c.liquidViscosity;
end
% Vectorized internal faces; identical harmonic conductances and upwinding.
k=(1:N-1)'; z=k+1;
G=conductance(p.k(k),p.k(k+1),g.dx(k),g.dx(k+1));
f.T(z)=G.*(s.T(k)-s.T(k+1));
ka=k(g.anode(k)&g.anode(k+1));
G=conductance(p.DH(ka),p.DH(ka+1),g.dx(ka),g.dx(ka+1));
f.H(ka+1)=G.*(p.cH(ka)-p.cH(ka+1));
kc=k(g.cathode(k)&g.cathode(k+1));
G=conductance(p.DO(kc),p.DO(kc+1),g.dx(kc),g.dx(kc+1));
f.O(kc+1)=G.*(p.cO(kc)-p.cO(kc+1));
kp=k(g.porous(k)&g.porous(k+1)&g.layer(k)~=4);
G=conductance(p.Dv(kp),p.Dv(kp+1),g.dx(kp),g.dx(kp+1));
f.v(kp+1)=G.*(p.rv(kp)-p.rv(kp+1));
if c.liquidMigration
    Gabs=conductance(K0(kp),K0(kp+1),g.dx(kp),g.dx(kp+1));
    dpc=pc(kp+1)-pc(kp); sUp=sl(kp);
    neg=dpc<0; sUp(neg)=sl(kp(neg)+1);
    f.l(kp+1)=c.rhoL/c.liquidViscosity*Gabs.*sUp.^3.*dpc;
end
ki=k(g.ion(k)&g.ion(k+1));
aL=g.omega(ki).^1.5*c.Cm.*max(p.Db(ki),0);
aR=g.omega(ki+1).^1.5*c.Cm.*max(p.Db(ki+1),0);
G=conductance(aL,aR,g.dx(ki),g.dx(ki+1));
lamUp=p.lambda(ki); if j<0,lamUp=p.lambda(ki+1);end
drag=2.5*max(lamUp,0)/22*c.Mw*j.*gammaFace(ki+1)/c.F;
f.b(ki+1)=G.*(p.lambda(ki)-p.lambda(ki+1))+drag;
Gt=1/(g.dx(1)/(2*p.k(1))+1/c.h);
f.T(1)=-Gt*(s.T(1)-d.Tamb);
Gt=1/(g.dx(end)/(2*p.k(end))+1/c.h);
f.T(end)=Gt*(s.T(end)-d.Tamb);
ia=find(g.anode,1,'first'); ic=find(g.cathode,1,'last');
f.H(ia)=2*p.DH(ia)/g.dx(ia)*(c.p/(c.R*s.T(ia))-p.cH(ia));
f.O(ic+1)=2*p.DO(ic)/g.dx(ic)*(p.cO(ic)-c.yO*c.p/(c.R*s.T(ic)));
f.v(ia)=-vaporOutlet(2*p.Dv(ia)/g.dx(ia)*p.rv(ia),j,s.T(ia),c,'H');
f.v(ic+1)=vaporOutlet(2*p.Dv(ic)/g.dx(ic)*p.rv(ic),j,s.T(ic),c,'O');
if c.liquidMigration
    f.l(ia)=2*mob(ia)/g.dx(ia)*pc(ia);
    f.l(ic+1)=-2*mob(ic)/g.dx(ic)*pc(ic);
end
end

function G=conductance(a,b,da,db)
G=zeros(size(a)); m=a>0&b>0;
G(m)=1./(da(m)./(2*a(m))+db(m)./(2*b(m)));
end

function J=vaporOutlet(Jdiff,j,T,c,side)
switch c.vaporBoundary
    case 'dry', J=Jdiff;
    case 'closed', J=0;
    case 'finite_flow'
        if side=='H', nDry=c.stoich*max(j,0)/(2*c.F);
        else, nDry=c.stoich*max(j,0)/(4*c.F*c.yO);end
        Tc=T-273.15;
        ps=611.15*exp((23.036-Tc/333.7)*Tc/(279.82+Tc));
        Jflow=nDry*c.Mw*ps/max(c.p-ps,1);
        J=max(Jdiff,0)*Jflow/(max(Jdiff,0)+Jflow+realmin);
    otherwise, error('Unknown vapor boundary');
end
end
