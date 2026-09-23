function c = q1_config()
% All lengths, currents and inventories use SI units.
c.F=96485.33212; c.R=8.314462618; c.Mw=0.018;
c.p=101325; c.yO=(0.233/31.998)/(0.233/31.998+0.767/28.014);
c.rhoL=990; c.rhoI=920; c.Cm=2150*c.Mw/1;
c.DHref=1.10e-4; c.DOref=3.20e-5; c.Dvref=2.82e-5;
c.kGas=.0246; c.kLiquid=.6; c.kIce=2.3;
c.omega=0.3; c.h=40; c.Lv=2.5e6; c.Lf=333600;
c.dxLayer=[.002 .00015 3.4e-6 12e-6 11.3e-6 .00015 .002];
c.CLayer=[1980*766 185*545 970*240 2150*1050 970*240 185*545 1980*766];
c.kLayer=[95 .3 .27 .24 .27 .3 95];
c.epsLayer=[0 .8 .3916 0 .4207 .8 0];
c.inputArea_cm2=303; c.inputMode='common_I_over_A';
c.j0ref=.01; c.Ea=67000; c.alpha=.5; c.fCL=1; c.Rcontact=1e-6;
c.rateBasis='bulk_Cm'; c.kdes=.001; c.kabs=1; c.kbi=1;
c.kvi=1e-4; c.kcond=1; c.kevap=1; c.kbl=.5;
c.kfreeze=1; c.kmelt=1; c.knf=1; c.kfn=1; c.donorM=1e-8;
c.vaporBoundary='finite_flow'; c.stoich=2;
c.concentrationModel='limiting_current'; c.concentrationFactor=1;
c.voltageResistance='power_average'; c.activitySurface='ice_below_zero';
c.reactionClosure='parallel_active'; % local_average retains v0.3 control
c.kineticHydrationExponent=0; % hypothesis-only sensitivity, not in baseline
c.temperatureAverage='seven_layer'; c.conductivity='springer1268';
c.iceEnabled=true; c.iceFeedback=true; c.pemFreeze=true;
c.phaseHeat=true; c.freezeConductivityLambda=false;
c.clFrozenWater=true; % v0.4 structural repair; false reproduces 272-state v0.3
c.icePartition='volume_weighted'; % alternative: pore_first
c.liquidRoute=false; c.liquidMigration=false;
c.freezingClosure='direct_ice'; % alternative: delayed_nucleation
c.nucleationThreshold=.2; c.nucleationRate=1;
c.surfaceTension=.075; c.liquidViscosity=.001;
c.K0Layer=[0 6.2e-12 6.2e-13 0 6.2e-13 6.2e-12 0];
c.thetaLayer=[0 110 100 0 100 110 0];
c.reltol=1e-6; c.abstol=1e-9; c.maxstep=.05;
c.runCalibration=false; c.version='q1-v0.4-implementation';
end
