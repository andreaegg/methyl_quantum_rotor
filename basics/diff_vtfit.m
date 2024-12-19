function diff = diff_vtfit(V3,vt,K,B)

vtcurr = hindered_rotor_model(V3,K,B);
vtcurr = vtcurr*boltzm/planck;
diff   = abs(vtcurr - vt);
end