function diff = diff_vtfit(V3,vt,m,n,B)

vtcurr = hindered_rotor_model(V3,m,n,B);
vtcurr = vtcurr*boltzm/planck;
diff   = abs(vtcurr - vt);
end