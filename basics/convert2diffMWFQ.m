function Bout = convert2diffMWFQ(mwfqin,Bin,mwfqout)
    
    g    = (planck*mwfqin*1e9)/(bmagn*Bin*1e-4);
    Bout = (planck*mwfqout*1e9)/(bmagn*g*1e-4);

end