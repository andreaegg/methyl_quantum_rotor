function [h,g,vt,V3] = fit_tunneldistribution(V,V0,sd,factor,nreigval,rotor)

vt = V3pot2vt(V,nreigval,rotor);
V3 = V;
g  = distribution_V3(V,V0,sd,factor);
h  = distribution_vt(g,vt.MHz,V);
end