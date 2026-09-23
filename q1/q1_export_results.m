function q1_export_results(out,folder,label)
if nargin<3,label='unfitted_preview';end
c=out.cfg; d=out.cases; runs=out.runs;
fid=fopen(fullfile(folder,'config.json'),'w');
fwrite(fid,jsonencode(c,'PrettyPrint',true)); fclose(fid);
names={'Faraday','gas_constant','water_molar_mass','ionomer_water_capacity', ...
    'pressure','oxygen_mole_fraction','convective_h','vapor_latent_heat', ...
    'fusion_latent_heat','input_area_m2','j0_ref','activation_energy','CL_factor', ...
    'H2_diffusivity_ref','O2_diffusivity_ref','vapor_diffusivity_ref', ...
    'liquid_density','ice_density'};
vals=[c.F;c.R;c.Mw;c.Cm;c.p;c.yO;c.h;c.Lv;c.Lf; ...
    c.inputArea_cm2*1e-4;c.j0ref;c.Ea;c.fCL;c.DHref;c.DOref;c.Dvref; ...
    c.rhoL;c.rhoI];
units={'C/mol';'J/(mol K)';'kg/mol';'kg/m3';'Pa';'1';'W/(m2 K)'; ...
    'J/kg';'J/kg';'m2';'A/m2';'J/mol';'1';'m2/s';'m2/s';'m2/s'; ...
    'kg/m3';'kg/m3'};
source={'physical_constant';'physical_constant';'rounded_attachment_convention'; ...
    'attachment_1';'attachment_1';'mass_to_mole_conversion';'attachment_1'; ...
    'attachment_1';'attachment_1';'review_input_reconstruction'; ...
    'attachment_1_default';'problem_default';'structural_parameter'; ...
    'working_assumption';'working_assumption';'working_assumption'; ...
    'attachment_1';'attachment_1'};
P=table(string(names(:)),vals,string(units(:)),string(source(:)), ...
    'VariableNames',{'parameter','value_SI','unit','source'});
writetable(P,fullfile(folder,'parameters_SI.csv'));
f=figure('Visible','off'); tiledlayout(3,1);
nexttile; hold on; for k=1:2,plot(d(k).t,d(k).j_Acm2);end
ylabel('j (A/cm^2)'); legend({d.name},'Location','best');
nexttile; hold on; for k=1:2,plot(d(k).t,d(k).V_exp);end
ylabel('V (V)');
nexttile; hold on; for k=1:2,plot(d(k).t,d(k).T_exp_C);end
ylabel('T (°C)'); xlabel('t (s)');
exportgraphics(f,fullfile(folder,'input_curves.png'),'Resolution',180); close(f);
f=figure('Visible','off'); tiledlayout(2,1);
nexttile; hold on;
for k=1:2
    plot(d(k).t,d(k).V_exp,'--');
    plot([runs(k).obs.time_s],[runs(k).obs.V_sim],'-');
end
ylabel('V (V)'); title(strrep(label,'_',' '));
nexttile; hold on;
for k=1:2
    plot(d(k).t,d(k).T_exp_C,'--');
    plot([runs(k).obs.time_s],[runs(k).obs.T_seven_C],'-');
end
ylabel('T (°C)'); xlabel('t (s)');
exportgraphics(f,fullfile(folder,'trajectory_preview.png'),'Resolution',180); close(f);
rows=[];
for k=1:2
    if ~strcmp(runs(k).status,'VALID'),continue;end
    obs=runs(k).obs; n=numel(obs);
    V=[obs.V_sim]'; TT=[obs.T_seven_C]';
    rows=[rows; table(string(d(k).name),string(label), ...
        sqrt(mean((V-d(k).V_exp).^2)),sqrt(mean((TT-d(k).T_exp_C).^2)), ...
        runs(k).balance.ew,runs(k).balance.eE, ...
        'VariableNames',{'case_name','result_type','V_RMSE','T_RMSE_K','water_error','energy_error'})]; %#ok<AGROW>
end
if ~isempty(rows),writetable(rows,fullfile(folder,'metrics.csv'));end
end
