function report=q1_local_phase_checks(root)
c=q1_config();c.freezingClosure='delayed_nucleation';c.nucleationMode='local';
c.nucleationThreshold=.03;c.liquidMigration=true;
[d,~]=q1_read_data(root,c);g=q1_make_grid(c,'coarse');ix=q1_make_index(g);
y=q1_initial(d(2),c,g,ix);target=find(g.cCL,1);
liquidIndex=find(ix.mlcells==target);y(ix.ml(liquidIndex))=.05*c.rhoL*g.eps(target);
[dy,diag]=q1_rhs_diagnostic(0,y,d(2),c,g,ix);
localIndex=find(ix.nuccells==target);
others=setdiff(1:numel(ix.nuc),localIndex);
rows={'cold_wet_cell_nucleates',dy(ix.nuc(localIndex))>0; ...
    'dry_remote_cells_do_not_nucleate',all(dy(ix.nuc(others))==0); ...
    'no_instant_ice_without_nuclei',all(diag.r.li==0)};
y(ix.T(target))=274;
[dy,~]=q1_rhs_diagnostic(0,y,d(2),c,g,ix);
rows(end+1,:)={'warm_wet_cell_does_not_nucleate',dy(ix.nuc(localIndex))==0};
% Add nuclei to the cold wet cell: freezing removes liquid and adds ice.
y(ix.T(target))=248.15;y(ix.nuc(localIndex))=.5;
[~,diag]=q1_rhs_diagnostic(0,y,d(2),c,g,ix);
rows(end+1,:)={'existing_local_nuclei_freeze_liquid',diag.r.li(target)>0};
report=cell2table(rows,'VariableNames',{'check','pass'});
writetable(report,fullfile(root,'results','q1_v05','local_phase_checks.csv'));
assert(all(report.pass),'Local phase check failed');
end
