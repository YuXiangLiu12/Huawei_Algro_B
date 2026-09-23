function g=q1_make_grid(c,level)
switch level
    case 'coarse', n=[2 6 3 6 6 6 2];
    case 'medium', n=[4 12 6 12 12 12 4];
    case 'fine', n=[8 24 12 24 24 24 8];
    case 'ion_refined', n=[8 24 24 48 48 24 8];
    otherwise, error('Unknown grid level');
end
g.level=level; g.layer=repelem((1:7)',n(:));
g.dx=repelem((c.dxLayer(:)./n(:)),n(:)); g.N=numel(g.dx);
g.xf=[0;cumsum(g.dx)]; g.xc=(g.xf(1:end-1)+g.xf(2:end))/2;
g.porous=ismember(g.layer,[2 3 5 6]);
g.ion=ismember(g.layer,[3 4 5]); g.pem=g.layer==4;
g.anode=ismember(g.layer,[2 3]); g.cathode=ismember(g.layer,[5 6]);
g.aCL=g.layer==3; g.cCL=g.layer==5; g.MEA=ismember(g.layer,[2 3 4 5 6]);
g.eps=c.epsLayer(g.layer)'; g.C=c.CLayer(g.layer)'; g.k=c.kLayer(g.layer)';
g.omega=zeros(g.N,1); g.omega(g.ion)=c.omega; g.omega(g.pem)=1;
g.LMEA=sum(g.dx(g.MEA)); g.LaCL=sum(g.dx(g.aCL)); g.LcCL=sum(g.dx(g.cCL));
g.gammaFace=zeros(g.N+1,1);
ia=find(g.aCL); ip=find(g.pem); ic=find(g.cCL);
g.gammaFace(ia(1):ia(end)+1)=linspace(0,1,numel(ia)+1);
g.gammaFace(ip(1):ip(end)+1)=1;
g.gammaFace(ic(1):ic(end)+1)=linspace(1,0,numel(ic)+1);
g.clFrozenWater=c.clFrozenWater;
g.freezingClosure=c.freezingClosure;
g.nucleationMode='global_legacy';
if isfield(c,'nucleationMode'),g.nucleationMode=c.nucleationMode;end
end
