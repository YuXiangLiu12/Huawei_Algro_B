function y=q1_initial(d,c,g,ix)
y=zeros(ix.n,1); y(ix.T)=d.T0;
y(ix.qH)=g.eps(g.anode)*c.p/(c.R*d.T0);
y(ix.qO)=g.eps(g.cathode)*c.yO*c.p/(c.R*d.T0);
y(ix.bn)=g.omega(g.ion)*c.Cm*3;
end
