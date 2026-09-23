function q1_plot_mesh(cases,runs,folder)
% Reusable plotter; runs may contain only the saved observation structures.
f=figure('Visible','off','Color','w','Position',[100 100 1200 900]);
cleanup=onCleanup(@()close(f));
tiledlayout(3,2,'TileSpacing','compact');
labels={'coarse','medium','fine'};colors=lines(3);
for k=1:numel(cases)
    d=cases(k);
    nexttile(k);hold on;
    plot(d.t,d.V_exp,'k.','DisplayName','Experiment');
    for z=1:3,plot(d.t,[runs{z}(k).obs.V_sim],'Color',colors(z,:),'DisplayName',labels{z});end
    title(sprintf('%g deg C',d.T0-273.15));ylabel('Voltage (V)');
    legend('Location','northeast');grid on;
    nexttile(k+2);hold on;
    plot(d.t,d.T_exp_C,'k.');
    for z=1:3,plot(d.t,[runs{z}(k).obs.T_seven_C],'Color',colors(z,:));end
    ylabel('Mean temperature (deg C)');grid on;
    nexttile(k+4);hold on;
    for z=1:3,plot(d.t,[runs{z}(k).obs.ice_volume_max],'Color',colors(z,:));end
    ylabel('Max. ice volume fraction (predicted)');xlabel('Time (s)');grid on;
end
exportgraphics(f,fullfile(folder,'mesh_comparison.png'),'Resolution',160);
end
