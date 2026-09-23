function o=q1_observe(t,y,d,c,g,ix)
s=q1_unpack(y,g,ix); p=q1_properties(s,c,g);
v=q1_voltage(s,p,d.jfun(t),c,g);
o.time_s=t; o.V_sim=v.V;
o.T_seven_C=sum(s.T.*g.dx)/sum(g.dx)-273.15;
o.T_MEA_C=sum(s.T(g.MEA).*g.dx(g.MEA))/g.LMEA-273.15;
[o.ice_volume_max,imax]=max(s.mi(g.porous)/c.rhoI);
porousCells=find(g.porous);o.ice_peak_layer=NaN;o.ice_peak_x_um=NaN;
if o.ice_volume_max>0
    o.ice_peak_layer=g.layer(porousCells(imax));
    o.ice_peak_x_um=1e6*g.xc(porousCells(imax));
end
o.ice_saturation_max=max(p.si(g.porous));
o.liquid_saturation_max=max(p.sl(g.porous));
o.lambda_aCL=mean(p.lambda(g.aCL));
o.lambda_PEM=mean(p.lambda(g.pem));
o.lambda_cCL=mean(p.lambda(g.cCL));
o.lambda_oversat_CL=max(p.lambda(g.aCL|g.cCL)-p.lambdaSat(g.aCL|g.cCL));
o.CL_frozen_ion_water=sum(s.bfc.*g.dx);
o.gas_porosity_min=min(p.eg(g.porous));
o.V_Erev=v.Erev; o.V_act=v.etaAct; o.V_ohm=v.etaOhm; o.V_con=v.etaCon;
o.R_aCL=v.RaCL; o.R_PEM=v.Rmem; o.R_cCL=v.RcCL;
o.jlim_H=v.jlimH; o.jlim_O=v.jlimO;
o.reaction_weight_max=max(v.reactionWeight);
o.nucleation_progress=max(s.nuc(:));
o.water_bound_kgm2=sum(s.bn.*g.dx);
o.water_vapor_kgm2=sum(s.mv.*g.dx);
o.water_liquid_kgm2=sum(s.ml.*g.dx);
o.water_pore_ice_kgm2=sum(s.mi.*g.dx);
o.water_ion_ice_kgm2=sum((s.bf+s.bfc).*g.dx);
end
