function v=q1_voltage(s,p,j,c,g)
Tc=sum(s.T(g.cCL).*g.dx(g.cCL))/g.LcCL;
v.Erev=1.229-8.5e-4*(Tc-298.15)+c.R*Tc/(2*c.F)* ...
    log((c.p/101325)*sqrt(c.yO*c.p/101325));
chi=max(1-p.si(g.cCL),0).^3.5;
if ~c.iceFeedback,chi=ones(size(chi));end
j0=c.j0ref*exp(-c.Ea/c.R*(1./s.T(g.cCL)-1/298.15));
j0=j0.*(max(p.lambda(g.cCL),.65)/3).^c.kineticHydrationExponent;
switch c.reactionClosure
    case 'local_average'
        v.etaAct=sum(g.dx(g.cCL).*c.R.*s.T(g.cCL)/(c.alpha*c.F).* ...
            asinh(j./(2*j0.*max(chi,1e-12))))/g.LcCL;
        v.reactionWeight=ones(sum(g.cCL),1);
        v.gammaFace=g.gammaFace;
    case 'parallel_active'
        j0eff=sum(g.dx(g.cCL).*j0.*chi)/g.LcCL;
        v.etaAct=c.R*Tc/(c.alpha*c.F)*asinh(j/(2*max(j0eff,1e-12)));
        if abs(j)<eps
            v.reactionWeight=ones(sum(g.cCL),1);
        else
            local=j0.*chi.*sinh(c.alpha*c.F*v.etaAct./(c.R*s.T(g.cCL)));
            meanLocal=sum(g.dx(g.cCL).*local)/g.LcCL;
            v.reactionWeight=local/max(meanLocal,realmin);
        end
        v.gammaFace=g.gammaFace;
        cells=find(g.cCL);
        v.gammaFace(cells(1))=1;
        for kk=1:numel(cells)
            cellID=cells(kk);
            v.gammaFace(cellID+1)=v.gammaFace(cellID)- ...
                v.reactionWeight(kk)*g.dx(cellID)/g.LcCL;
        end
    otherwise,error('Unknown reaction closure');
end
kappa=p.kappa;
if c.freezeConductivityLambda
    lam=3; T=s.T(g.ion);
    if strcmp(c.conductivity,'cold2222')
        kap=100*exp(2222*(1/303-1./T))*(.005139*lam-.00326);
    else
        kap=(.5139*lam-.326).*exp(1268*(1/303.15-1./T));
    end
    kap=kap.*(1+(c.fCL*g.omega(g.ion).^1.5-1).*(g.layer(g.ion)~=4));
    kappa(g.ion)=kap;
end
gam=v.gammaFace;
switch c.voltageResistance
    case 'power_average'
        weight=(gam(1:end-1).^2+gam(1:end-1).*gam(2:end)+gam(2:end).^2)/3;
    case 'endpoint_both_CL'
        weight=(gam(1:end-1)+gam(2:end))/2;
    case 'mao_cathode_only'
        weight=(gam(1:end-1)+gam(2:end))/2;
        weight(g.aCL)=0;
    otherwise,error('Unknown voltage resistance closure');
end
Rk=zeros(g.N,1); Rk(g.ion)=g.dx(g.ion).*weight(g.ion)./max(kappa(g.ion),1e-12);
v.RaCL=sum(Rk(g.aCL)); v.Rmem=sum(Rk(g.pem)); v.RcCL=sum(Rk(g.cCL));
v.etaOhm=j*(sum(Rk)+c.Rcontact);
% Problem Eq. (43)-(44): oxygen limiting current based on reaction-zone
% concentration and effective diffusion distance. Series resistance gives
% D_eff/L_diff without inventing one layer's diffusivity as the whole path.
ih=find(g.anode); io=find(g.cathode);
RH=sum(g.dx(ih)./max(p.DH(ih),1e-20));
RO=sum(g.dx(io)./max(p.DO(io),1e-20));
ch=sum(p.cH(g.aCL).*g.dx(g.aCL))/g.LaCL;
co=sum(p.cO(g.cCL).*g.dx(g.cCL))/g.LcCL;
v.jlimH=2*c.F*ch/RH;
v.jlimO=4*c.F*co/RO;
if strcmp(c.concentrationModel,'limiting_current')
    v.etaCon=-c.concentrationFactor*c.R*Tc/(4*c.F)* ...
        log(max(1-j/max(v.jlimO,1e-12),1e-12));
elseif strcmp(c.concentrationModel,'reaction_concentration')
    ph=sum(p.cH(g.aCL).*s.T(g.aCL).*g.dx(g.aCL))/g.LaCL*c.R;
    po=sum(p.cO(g.cCL).*s.T(g.cCL).*g.dx(g.cCL))/g.LcCL*c.R;
    v.etaCon=c.R*Tc/(2*c.F)*log(c.p/max(ph,1e-12)) ...
        +c.R*Tc/(4*c.F)*log(c.yO*c.p/max(po,1e-12));
else,error('Unknown concentration model');end
v.V=v.Erev-v.etaAct-v.etaOhm-v.etaCon;
end
