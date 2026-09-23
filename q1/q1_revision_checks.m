function report=q1_revision_checks(root)
% Regression test of face vectorization against the preserved v0.4 code.
if nargin<1,root=fileparts(fileparts(mfilename('fullpath')));end
folder=fullfile(root,'results','q1_v05');checkdir=fullfile(folder,'checkcode');
if ~isfolder(checkdir),mkdir(checkdir);end
src=fileread(fullfile(folder,'source_before','q1_face_fluxes.m'));
src=strrep(src,'function f=q1_face_fluxes(','function f=q1_face_fluxes_v04_reference(');
fid=fopen(fullfile(checkdir,'q1_face_fluxes_v04_reference.m'),'w');fwrite(fid,src);fclose(fid);
addpath(checkdir);cleanup=onCleanup(@()rmpath(checkdir));
c=q1_config();[cases,~]=q1_read_data(root,c);g=q1_make_grid(c,'coarse');ix=q1_make_index(g);
y=q1_initial(cases(2),c,g,ix);s=q1_unpack(y,g,ix);
s.T=s.T+linspace(0,4,g.N)';
s.ml(g.porous)=c.rhoL*g.eps(g.porous).*linspace(.01,.12,sum(g.porous))';
s.mi(g.porous)=c.rhoI*g.eps(g.porous).*linspace(.03,.08,sum(g.porous))';
s.mv(g.porous)=linspace(1e-5,1e-4,sum(g.porous))';
rows=cell(0,4);
for liquid=[false true]
    cc=c;cc.liquidMigration=liquid;p=q1_properties(s,cc,g);
    for j=[0 10 2000]
        v=q1_voltage(s,p,j,cc,g);
        a=q1_face_fluxes(s,p,j,cases(2),cc,g,v);
        b=q1_face_fluxes_v04_reference(s,p,j,cases(2),cc,g,v);
        names=fieldnames(a);err=0;
        for z=1:numel(names)
            err=max(err,max(abs(a.(names{z})-b.(names{z}))./max(1,abs(b.(names{z})))));
        end
        rows(end+1,:)={liquid,j,err,err<1e-12};
    end
end
report=cell2table(rows,'VariableNames',{'liquid_migration','current_SI','relative_error','pass'});
writetable(report,fullfile(folder,'flux_regression.csv'));
assert(all(report.pass),'Vectorized flux regression failed');
end
