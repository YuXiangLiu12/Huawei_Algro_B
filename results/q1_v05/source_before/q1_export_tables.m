function q1_export_tables(out,folder,resultType)
% Call only for a declared calibrated fit and complete valid trajectories.
assert(all(strcmp({out.runs.status},'VALID')),'Full valid trajectories required');
assert(~strcmp(resultType,'unfitted_preview'),'Result type must disclose calibration');
times=(0:5:35)';
for k=1:2
    d=out.cases(k); run=out.runs(k); ix=run.index; g=run.grid;
    obs=repmat(run.obs(1),numel(times),1); Ve=zeros(size(times)); Te=Ve;
    for i=1:numel(times)
        y=deval(run.sol,times(i)).*ix.scale;
        obs(i)=q1_observe(times(i),y,d,out.cfg,g,ix);
        [delta,z]=min(abs(d.t-times(i)));
        assert(delta<1e-9,'Experimental time not found within tolerance');
        Ve(i)=d.V_exp(z); Te(i)=d.T_exp_C(z);
    end
    V=[obs.V_sim]'; T=[obs.T_seven_C]'; ice=[obs.ice_volume_max]';
    tab=table(times,V,Ve,100*abs(V-Ve)./abs(Ve),T,Te, ...
        100*abs(T-Te)./abs(Te),ice,[obs.ice_saturation_max]', ...
        repmat(string(resultType),numel(times),1), ...
        'VariableNames',{'time_s','V_sim','V_exp','V_error_pct', ...
        'T_sim_C','T_exp_C','T_error_pct','ice_absolute_max', ...
        'ice_pore_saturation_max','result_type'});
    writetable(tab,fullfile(folder,sprintf('table%d.csv',k)));
end
end
