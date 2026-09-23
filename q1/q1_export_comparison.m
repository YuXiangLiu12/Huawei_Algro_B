function q1_export_comparison(root,out,folder)
% Direct source-data comparison; ice curves are predictions, not observations.
cols=[.10 .36 .70;.78 .25 .14];
f=figure('Visible','off','Color','w','Position',[100 100 1180 930]);
tl=tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
for k=1:2
    d=out.cases(k);o=struct2table(out.runs(k).obs);
    old=readtable(fullfile(root,'results','q1_v04',sprintf('case%d_observations.csv',k)));
    nexttile(k);hold on;
    plot(d.t,d.V_exp,'k.','MarkerSize',6,'DisplayName','Experiment');
    plot(old.time_s,old.V_sim,'--','Color',[.6 .6 .6],'LineWidth',1.2,'DisplayName','Original');
    plot(o.time_s,o.V_sim,'Color',cols(k,:),'LineWidth',1.8,'DisplayName','Calibrated');
    ylabel('Voltage (V)');title(sprintf('%d deg C',-15-5*k));
    legend('Location','best');box off;grid on;xlim([0 36.6]);
    nexttile(k+2);hold on;
    plot(d.t,d.T_exp_C,'k.','MarkerSize',6);
    plot(old.time_s,old.T_seven_C,'--','Color',[.6 .6 .6],'LineWidth',1.2);
    plot(o.time_s,o.T_seven_C,'Color',cols(k,:),'LineWidth',1.8);
    ylabel('Mean temperature (deg C)');box off;grid on;xlim([0 36.6]);
    nexttile(k+4);hold on;
    plot(old.time_s,old.ice_volume_max,'--','Color',[.6 .6 .6],'LineWidth',1.2);
    plot(o.time_s,o.ice_volume_max,'Color',cols(k,:),'LineWidth',1.8);
    ylabel('Maximum ice volume fraction');xlabel('Time (s)');
    box off;grid on;xlim([0 36.6]);ylim([0 max(.02,max(old.ice_volume_max)*1.08)]);
end
set(findall(f,'Type','axes'),'FontName','Arial','FontSize',11);
title(tl,'Q1: joint calibration to both temperatures (ice is unobserved)');
exportgraphics(f,fullfile(folder,'comparison.png'),'Resolution',180);
exportgraphics(f,fullfile(folder,'comparison.pdf'),'ContentType','vector');close(f);
end
