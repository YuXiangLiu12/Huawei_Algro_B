function run=q1_simulate(d,c,level)
if nargin<3,level='medium';end
g=q1_make_grid(c,level); ix=q1_make_index(g); y0=q1_initial(d,c,g,ix);
scale=ix.scale(:); z0=y0./scale;
rhs=@(t,z) q1_rhs(t,z.*scale,d,c,g,ix)./scale;
ev=@(t,z) q1_events(t,z.*scale,d,c,g,ix);
opt=odeset('RelTol',c.reltol,'AbsTol',c.abstol,'MaxStep',c.maxstep, ...
    'NonNegative',ix.nonnegative,'Events',ev);
run.status='VALID'; run.event='';
try
    sol=ode15s(rhs,[d.t(1) d.t(end)],z0,opt);
    run.t_end=sol.x(end);
    if run.t_end<d.t(end)-1e-7
        run.status='EARLY_STOP';
        labels={'PORE_CLOSED','OXYGEN_DEPLETED','HYDROGEN_DEPLETED', ...
            'DRY_MODEL_LIMIT','DIFFUSIVITY_LIMIT','CONDUCTIVITY_LIMIT', ...
            'OXYGEN_LIMITING_CURRENT'};
        if isfield(sol,'ie') && ~isempty(sol.ie)
            run.event=labels{sol.ie(end)};
        else
            run.event='UNKNOWN_EVENT';
        end
    end
    ts=d.t(d.t<=run.t_end+1e-9); run.sample_t=ts;
    Y=deval(sol,ts).*scale;
    run.balance=q1_balance(Y,c,g,ix);
    if strcmp(run.status,'VALID') && (run.balance.ew>=1e-4 || run.balance.eE>=1e-3)
        run.status='BALANCE_FAIL';
    end
    run.Y=Y; run.sol=sol; run.grid=g; run.index=ix;
    run.obs=repmat(q1_observe(ts(1),Y(:,1),d,c,g,ix),numel(ts),1);
    for k=1:numel(ts)
        run.obs(k)=q1_observe(ts(k),Y(:,k),d,c,g,ix);
    end
catch ME
    run.status='SOLVER_ERROR'; run.event=ME.message; run.t_end=NaN;
    run.grid=g; run.index=ix;
end
end
