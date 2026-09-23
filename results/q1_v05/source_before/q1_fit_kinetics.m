function fit=q1_fit_kinetics(cases,c)
% Low-current 0--6 s approximation, joint two-temperature estimation.
% Unknowns are log(j0 at 253.15 K) and Ea (kJ/mol).
g=q1_make_grid(c,'medium'); ix=q1_make_index(g);
T=[]; J=[]; V=[];
for k=1:numel(cases)
    d=cases(k); m=d.t<=6+1e-9;
    T=[T; d.T_exp_K(m)]; J=[J; d.j_SI(m)]; V=[V; d.V_exp(m)];
end
theta0=[log(.01*exp(-c.Ea/c.R*(1/253.15-1/298.15))); 67];
obj=@(th) sum((predict(th,T,J,c,g,ix)-V).^2);
opt=optimset('Display','off','TolX',1e-8,'TolFun',1e-12);
th=fminsearch(@(z) obj([z(1);26+100/(1+exp(-z(2)))]), ...
    [theta0(1);log((67-26)/(126-67))],opt);
th=[th(1);26+100/(1+exp(-th(2)))];
fit.j0_253=exp(th(1)); fit.Ea_Jmol=th(2)*1000;
fit.j0_298=fit.j0_253*exp(fit.Ea_Jmol/c.R*(1/253.15-1/298.15));
fit.rmse_V=sqrt(obj(th)/numel(V));
fixed=fminsearch(@(z)obj([z;67]),theta0(1),opt);
fit.j0_253_Ea67=exp(fixed);
fit.rmse_Ea67_V=sqrt(obj([fixed;67])/numel(V));
fit.window_s=[0 6]; fit.note='Joint two-temperature low-current approximation; not a full-model fit';
end

function pred=predict(th,T,J,c,g,ix)
pred=zeros(size(T)); cc=c; cc.Ea=th(2)*1000;
cc.j0ref=exp(th(1))*exp(cc.Ea/cc.R*(1/253.15-1/298.15));
for i=1:numel(T)
    d.T0=T(i); d.Tamb=T(i);
    y=q1_initial(d,cc,g,ix); s=q1_unpack(y,g,ix);
    p=q1_properties(s,cc,g); v=q1_voltage(s,p,J(i),cc,g);
    pred(i)=v.V;
end
end
