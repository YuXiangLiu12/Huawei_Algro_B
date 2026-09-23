function checks=q1_checks(root)
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
c=q1_config(); [d,a]=q1_read_data(root,c);
names={}; values=[]; limits=[];
g=q1_make_grid(c,'medium'); ix=q1_make_index(g);
record('U01_length_m',abs(sum(g.dx)-.0043267),1e-12);
record('U01_heat_capacity',abs(sum(g.C.*g.dx)-6127.48),.05);
record('U01_rebuilt_charge',abs(a(1).charge_Ccm2-6.60221122),1e-7);
record('U01_272_state_compatibility',abs(compatStateCount(c)-272),0);
record('U01_290_state_repair',abs(ix.nphysical-290),0);
record('U01_shared_current',max(abs(d(1).j_SI-d(2).j_SI)),1e-8);
cc=c; cc.iceEnabled=false; cc.kdes=0; cc.kabs=0;
cc.vaporBoundary='closed'; cc.h=0;
g0=q1_make_grid(cc,'medium'); ix0=q1_make_index(g0);
dd=d(1); dd.jfun=@(t)0; y=q1_initial(dd,cc,g0,ix0);
dy=q1_rhs(0,y,dd,cc,g0,ix0);
record('U02_closed_uniform_stationary',max(abs(dy)),1e-6);
% Exact two-control-volume internal flux cancellation, independent of dx.
aa=.2; bb=.7; dxL=.001; dxR=.003; uL=5; uR=2;
J=(uL-uR)/(dxL/(2*aa)+dxR/(2*bb));
record('U04_face_conservation',abs((-J/dxL)*dxL+(J/dxR)*dxR),1e-12);
% Seven transfers conserve total water by column construction.
A=[-1 -1 -1 0 0 0 -1 -1;0 0 0 0 0 0 1 0; ...
    1 0 0 -1 -1 0 0 0;0 1 0 1 0 -1 0 0; ...
    0 0 1 0 1 1 0 0;0 0 0 0 0 0 0 1];
record('U06_water_transfer_columns',max(abs(sum(A,1))),0);
% Faraday source, integrated independently over each CL.
j=2000; WH=c.Mw*j/(2*c.F*g.LcCL)*sum(g.dx(g.cCL));
record('U05_faraday_water',abs(WH-c.Mw*j/(2*c.F)),1e-13);
% Frozen/melted mass and latent heat are exact scalar identities.
m0=.05; t=.7; mf=m0*exp(-c.kfreeze*t);
record('U07_freezing_analytic',abs((m0-mf)-(m0*(1-exp(-c.kfreeze*t)))),1e-14);
record('U08_latent_sign',double(~(c.Lf*c.kfreeze*m0>0 && ...
    c.Lf*(-c.kmelt*m0)<0)),0);
% The membrane flux has no boundary exit, including drag at j>0.
s=q1_unpack(y,g0,ix0); p=q1_properties(s,cc,g0);
f=q1_face_fluxes(s,p,j,dd,cc,g0);
div=(f.b(1:end-1)-f.b(2:end))./g0.dx;
record('U09_ion_flux_conservation',abs(sum(div(g0.ion).*g0.dx(g0.ion))),1e-12);
% For uniform conductivity the exact CL power-average is L/(3*kappa).
kap=p.kappa(g0.cCL); gamma2=(g0.gammaFace(1:end-1).^2+ ...
    g0.gammaFace(1:end-1).*g0.gammaFace(2:end)+g0.gammaFace(2:end).^2)/3;
Rc=sum(g0.dx(g0.cCL).*gamma2(g0.cCL)./kap);
record('U11_cCL_resistance',abs(Rc-g0.LcCL/(3*kap(1))),1e-12);
% Ice heterogeneity must redistribute, never change total imposed current.
cc.reactionClosure='parallel_active';
s.mi(g0.cCL)=linspace(0,.2,numel(find(g0.cCL)))'.*cc.rhoI.*g0.eps(g0.cCL);
p=q1_properties(s,cc,g0); v=q1_voltage(s,p,j,cc,g0);
record('U13_reaction_weight_conservation', ...
    abs(sum(v.reactionWeight.*g0.dx(g0.cCL))/g0.LcCL-1),1e-12);
cells=find(g0.cCL);
record('U13_proton_current_boundaries', ...
    abs(v.gammaFace(cells(1))-1)+abs(v.gammaFace(cells(end)+1)),1e-12);
% Wet-to-dry liquid transport must not be blocked by dry receiver mobility.
cc.liquidMigration=true; cells=find(g0.layer==6); k=cells(2);
s.ml(k)=.2*cc.rhoL*g0.eps(k);
s.ml(k+1)=.01*cc.rhoL*g0.eps(k+1);
p=q1_properties(s,cc,g0); f=q1_face_fluxes(s,p,j,dd,cc,g0);
record('U14_liquid_flux_direction',double(f.l(k+1)<=0),0);
checks=table(string(names),values,limits,values<=limits, ...
    'VariableNames',{'check','error','limit','pass'});
folder=fullfile(root,'results','checks'); if ~isfolder(folder),mkdir(folder);end
writetable(checks,fullfile(folder,'q1_checks.csv'));
assert(all(checks.pass),'One or more module checks failed');
    function record(name,value,limit)
        names{end+1,1}=name; values(end+1,1)=value; limits(end+1,1)=limit;
    end
end

function n=compatStateCount(c)
c.clFrozenWater=false; g=q1_make_grid(c,'medium'); ix=q1_make_index(g); n=ix.nphysical;
end
