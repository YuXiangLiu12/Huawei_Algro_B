function r=q1_phase_rates(s,p,c,g)
N=g.N; names={'bv','bl','bi','vl','vi','li','bf','bfc'};
for i=1:numel(names),r.(names{i})=zeros(N,1);end
m=g.aCL|g.cCL; cold=s.T<273.15;
gate=max(1-p.sl-p.si,0);
donorV=max(s.mv,0)./(max(s.mv,0)+c.donorM);
lambdaEq=p.lambdaEq;
r.bv(m)=c.Cm*gate(m).*(c.kdes*max(p.lambda(m)-lambdaEq(m),0) ...
    -c.kabs*max(lambdaEq(m)-p.lambda(m),0).*donorV(m));
if strcmp(c.freezingClosure,'delayed_nucleation')
    if ~c.iceEnabled,error('Delayed nucleation requires iceEnabled');end
    mc=m&cold;
    r.bl(mc)=c.kbl*c.Cm.*max(p.lambda(mc)-p.lambdaSat(mc),0).*gate(mc);
    mw=g.porous;
    donorL=max(s.ml,0)./(max(s.ml,0)+c.donorM);
    r.vl(mw)=max(p.eg(mw),0)*c.Mw./(c.R*s.T(mw)).* ...
        (c.kcond*max(p.pv(mw)-p.psl(mw),0)- ...
        c.kevap*max(p.psl(mw)-p.pv(mw),0).*donorL(mw));
    mp=g.porous&cold;
    r.li(mp)=c.kfreeze*max(s.nuc,0).*max(s.ml(mp),0);
    mp=g.porous&~cold;
    r.li(mp)=-c.kmelt*max(s.mi(mp),0);
    if c.pemFreeze
        mp=g.pem&cold;
        r.bf(mp)=max(s.nuc,0)*c.knf*c.Cm.* ...
            max(p.lambda(mp)-p.lambdaSat(mp),0);
        mp=g.pem&~cold; r.bf(mp)=-c.kfn*max(s.bf(mp),0);
    end
    return;
end
if c.iceEnabled
    mc=m&cold;
    excess=c.kbi*c.Cm.*max(p.lambda(mc)-p.lambdaSat(mc),0);
    if g.clFrozenWater
        if strcmp(c.icePartition,'volume_weighted')
            poreShare=max(p.eg(mc),0)./(max(p.eg(mc),0)+g.omega(mc)).*gate(mc).^2;
        elseif strcmp(c.icePartition,'pore_first')
            poreShare=gate(mc).^2;
        else
            error('Unknown CL ice partition');
        end
        r.bi(mc)=excess.*poreShare;
        r.bfc(mc)=excess.*(1-poreShare);
    else
        r.bi(mc)=excess.*gate(mc);
    end
    mp=g.porous&cold;
    r.vi(mp)=c.kvi*max(p.eg(mp),0)*c.Mw./(c.R*s.T(mp)).*max(p.pv(mp)-p.psi(mp),0);
end
if g.clFrozenWater
    mc=m&~cold; r.bfc(mc)=-c.kfn*max(s.bfc(mc),0);
end
mw=g.porous&(~cold|c.liquidRoute);
donorL=max(s.ml,0)./(max(s.ml,0)+c.donorM);
r.vl(mw)=max(p.eg(mw),0)*c.Mw./(c.R*s.T(mw)).* ...
    (c.kcond*max(p.pv(mw)-p.psl(mw),0)- ...
    c.kevap*max(p.psl(mw)-p.pv(mw),0).*donorL(mw));
mb=m&(~cold|c.liquidRoute)&p.activity>=1;
r.bl(mb)=c.kbl*c.Cm.*max(p.lambda(mb)-lambdaEq(mb),0).*gate(mb);
mp=g.porous;
r.li(mp&cold&c.iceEnabled)=c.kfreeze*max(s.ml(mp&cold&c.iceEnabled),0);
r.li(mp&~cold)=-c.kmelt*max(s.mi(mp&~cold),0);
if c.pemFreeze && c.iceEnabled
    m=g.pem&cold;
    r.bf(m)=c.knf*c.Cm.*max(p.lambda(m)-p.lambdaSat(m),0) ...
        -c.kfn*min(max(s.bf(m),0),c.Cm*max(p.lambdaSat(m)-p.lambda(m),0));
    m=g.pem&~cold; r.bf(m)=-c.kfn*max(s.bf(m),0);
end
end
