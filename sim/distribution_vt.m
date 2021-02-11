function h = distribution_vt(g,vt,V3)

dV3dvt = gradient(V3)./gradient(vt);
h = -g.*dV3dvt;
h = h/sum(h);
end