function [cases,audit] = q1_read_data(root,cfg)
% Read original observations; rebuild current density from the unsmoothed I column.
f=fullfile(root,'氢燃料电池低温冷启动建模与控制策略研究  附件','附件2.xlsx');
sh=sheetnames(f); assert(numel(sh)==2,'Expected two worksheets');
cases=repmat(struct(),1,2);
for k=1:2
    M=readmatrix(f,'Sheet',sh{k},'Range','A3:E186');
    assert(isequal(size(M),[184 5]) && all(isfinite(M(:))));
    t=M(:,1); assert(all(diff(t)>0) && abs(t(end)-36.6)<1e-8);
    cases(k).name=sh{k}; cases(k).t=t; cases(k).I_A=M(:,2);
    cases(k).V_exp=M(:,3); cases(k).T_exp_C=M(:,4);
    cases(k).T_exp_K=M(:,4)+273.15; cases(k).j_report_Acm2=M(:,5);
    cases(k).T0=258.15-5*k; cases(k).Tamb=cases(k).T0;
    cases(k).j_Acm2=cases(k).I_A/cfg.inputArea_cm2;
    cases(k).j_SI=1e4*cases(k).j_Acm2;
    cases(k).jfun=griddedInterpolant(t,cases(k).j_SI,'linear','none');
    Q=trapz(t,cases(k).j_Acm2);
    audit(k).name=sh{k}; audit(k).N=numel(t); audit(k).area_cm2=cfg.inputArea_cm2;
    audit(k).charge_Ccm2=Q; audit(k).reportedCharge_Ccm2=trapz(t,M(:,5));
    audit(k).water_mgcm2=1e6*cfg.Mw/(2*cfg.F)*Q;
    audit(k).j_report_max_relative_difference=max(abs(M(:,5)-cases(k).j_Acm2)./cases(k).j_Acm2);
    audit(k).Vmin=min(M(:,3)); audit(k).Tfinal_C=M(end,4);
end
assert(max(abs(cases(1).I_A-cases(2).I_A))<1e-10,'Current programs differ');
end
