function ix=q1_make_index(g)
n=0; names={'T','qH','qO','mv','ml','mi','bn','bf'};
masks={true(g.N,1),g.anode,g.cathode,g.porous,g.porous,g.porous,g.ion,g.pem};
scales=[300 50 50 .001 100 100 100 100];
for k=1:numel(names)
    z=find(masks{k}); ix.(names{k})=(n+1:n+numel(z))';
    ix.([names{k} 'cells'])=z; n=n+numel(z);
    ix.scale(ix.(names{k}),1)=scales(k);
end
ix.bfc=[]; ix.bfccells=[];
if g.clFrozenWater
    ix.bfccells=find(g.aCL|g.cCL);
    ix.bfc=(n+1:n+numel(ix.bfccells))'; n=n+numel(ix.bfccells);
    ix.scale(ix.bfc,1)=100;
end
ix.nuc=[]; ix.nuccells=[];
if strcmp(g.freezingClosure,'delayed_nucleation')
    if strcmp(g.nucleationMode,'local')
        ix.nuccells=find(g.porous);
        ix.nuc=(n+1:n+numel(ix.nuccells))';n=n+numel(ix.nuccells);
    else
        ix.nuc=n+1;n=n+1;
    end
    ix.scale(ix.nuc,1)=1;
end
ix.Q=n+1; ix.Mout=n+2; ix.Mvout=n+3;
ix.Eec=n+4; ix.Econv=n+5; ix.Ephase=n+6;
ix.scale(n+1:n+6,1)=[1e4;1e-3;1e-3;1e4;1e4;1e4];
ix.n=n+6; ix.nphysical=n;
ix.nonnegative=[ix.qH;ix.qO;ix.mv;ix.ml;ix.mi;ix.bn;ix.bf;ix.bfc;ix.nuc;ix.Q];
end
