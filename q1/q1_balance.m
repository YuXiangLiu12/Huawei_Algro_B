function b=q1_balance(Y,c,g,ix)
n=size(Y,2); M=zeros(1,n); H=M; ice=M;
for k=1:n
    s=q1_unpack(Y(:,k),g,ix);
    M(k)=sum((s.bn+s.bf+s.bfc+s.mv+s.ml+s.mi).*g.dx);
    H(k)=sum(g.C.*(s.T-273.15).*g.dx)+c.Lv*sum(s.mv.*g.dx) ...
        -c.Lf*sum((s.mi+s.bf+s.bfc).*g.dx);
    ice(k)=max(s.mi(g.porous)/c.rhoI);
end
b.waterInventory=M; b.enthalpy=H; b.iceMax=ice;
b.waterResidual=M-M(1)-c.Mw/(2*c.F)*Y(ix.Q,:)+Y(ix.Mout,:);
b.energyResidual=H-H(1)-Y(ix.Eec,:)+Y(ix.Econv,:)+c.Lv*Y(ix.Mvout,:);
b.ew=max(abs(b.waterResidual))/max([M(1),c.Mw/(2*c.F)*Y(ix.Q,end),1e-6]);
b.eE=max(abs(b.energyResidual))/max(abs(Y(ix.Eec,end)),1e3);
end
