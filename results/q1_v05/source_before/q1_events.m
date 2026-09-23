function [value,isterminal,direction]=q1_events(t,y,d,c,g,ix)
s=q1_unpack(y,g,ix); p=q1_properties(s,c,g);
j=d.jfun(t); v=q1_voltage(s,p,j,c,g);
co=c.yO*c.p/(c.R*mean(s.T(g.cCL)));
ch=c.p/(c.R*mean(s.T(g.aCL)));
value=[min(p.eg(g.porous)./g.eps(g.porous))-1e-6; ...
    min(p.cO(g.cCL))/co-1e-6; min(p.cH(g.aCL))/ch-1e-6; ...
    min(p.lambda(g.ion))-.65; min(p.Db(g.ion))-1e-15; ...
    min(p.kappa(g.ion))-1e-8; v.jlimO/max(j,eps)-1];
isterminal=ones(numel(value),1); direction=-ones(numel(value),1);
end
