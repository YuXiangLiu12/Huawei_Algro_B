function s=q1_unpack(y,g,ix)
s.T=y(ix.T); vars={'qH','qO','mv','ml','mi','bn','bf','bfc'};
for k=1:numel(vars)
    name=vars{k}; s.(name)=zeros(g.N,1);
    tmp=s.(name); tmp(ix.([name 'cells']))=y(ix.(name)); s.(name)=tmp;
end
if isempty(ix.nuc)
    s.nuc=1;
elseif strcmp(g.nucleationMode,'local')
    s.nuc=zeros(g.N,1);s.nuc(ix.nuccells)=y(ix.nuc);
else
    s.nuc=y(ix.nuc);
end
end
