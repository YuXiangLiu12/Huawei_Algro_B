function [dy,diag]=q1_rhs_diagnostic(t,y,d,c,g,ix)
assert(numel(y)==ix.n && all(isfinite(y)),'Invalid state');
s=q1_unpack(y,g,ix); j=d.jfun(t);
p=q1_properties(s,c,g); v=q1_voltage(s,p,j,c,g);
f=q1_face_fluxes(s,p,j,d,c,g,v); r=q1_phase_rates(s,p,c,g);
% Rows: bound, PEM ice, vapor, liquid, pore ice. Each column sums to zero.
A=[-1 -1 -1 0 0 0 -1 -1;0 0 0 0 0 0 1 0; ...
    1 0 0 -1 -1 0 0 0;0 1 0 1 0 -1 0 0; ...
    0 0 1 0 1 1 0 0;0 0 0 0 0 0 0 1];
rate=[r.bv r.bl r.bi r.vl r.vi r.li r.bf r.bfc]';
S=A*rate; N=g.N;
div=@(J) (J(1:N)-J(2:N+1))./g.dx;
dt=div(f.T); phase=-c.Lv*r.bv+c.Lv*r.vl+(c.Lv+c.Lf)*r.vi ...
    +c.Lf*(r.bi+r.li+r.bf+r.bfc);
if ~c.phaseHeat,phase(:)=0;end
ec=zeros(N,1); ec(g.MEA)=j*(1.48-v.V)/g.LMEA;
dt=(dt+phase+ec)./g.C;
dH=div(f.H); dH(g.aCL)=dH(g.aCL)-j/(2*c.F*g.LaCL);
dO=div(f.O); dO(g.cCL)=dO(g.cCL)- ...
    j*v.reactionWeight/(4*c.F*g.LcCL);
db=div(f.b)+S(1,:)'; db(g.cCL)=db(g.cCL)+ ...
    c.Mw*j*v.reactionWeight/(2*c.F*g.LcCL);
dv=div(f.v)+S(3,:)'; dl=div(f.l)+S(4,:)'; di=S(5,:)'; df=S(2,:)'; dbc=S(6,:)';
dy=zeros(ix.n,1); dy(ix.T)=dt;
dy(ix.qH)=dH(g.anode); dy(ix.qO)=dO(g.cathode);
dy(ix.mv)=dv(g.porous); dy(ix.ml)=dl(g.porous); dy(ix.mi)=di(g.porous);
dy(ix.bn)=db(g.ion); dy(ix.bf)=df(g.pem);
if ~isempty(ix.bfc),dy(ix.bfc)=dbc(g.aCL|g.cCL);end
if ~isempty(ix.nuc)
    dy(ix.nuc)=c.nucleationRate*max(max(p.sl(g.porous))-c.nucleationThreshold,0) ...
        *max(1-s.nuc,0);
end
ia=find(g.anode,1,'first'); ic=find(g.cathode,1,'last');
dy(ix.Q)=j; dy(ix.Mout)=f.v(ic+1)-f.v(ia)+f.l(ic+1)-f.l(ia);
dy(ix.Mvout)=f.v(ic+1)-f.v(ia);
dy(ix.Eec)=j*(1.48-v.V);
dy(ix.Econv)=f.T(end)-f.T(1);
dy(ix.Ephase)=sum(phase.*g.dx);
diag.s=s; diag.p=p; diag.f=f; diag.r=r; diag.v=v;
end
