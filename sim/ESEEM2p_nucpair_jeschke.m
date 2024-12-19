function Vori = ESEEM2p_nucpair_jeschke(tau,Sys,ori)

    yI   = nucgval('1H')*nmagn/hbar;
    ye   = Sys.g*bmagn/hbar;
    rdd1 = norm(Sys.Scoord-Sys.Icoord(1,:));
    wdd1 = ye*yI*mu0*hbar/(4*pi*(rdd1*1e-10)^3)/(2*pi*1e6);
    uv1  = (Sys.Scoord-Sys.Icoord(1,:))/rdd1;
    rdd2 = norm(Sys.Scoord-Sys.Icoord(2,:));
    wdd2 = ye*yI*mu0*hbar/(4*pi*(rdd2*1e-10)^3)/(2*pi*1e6);
    uv2  = (Sys.Scoord-Sys.Icoord(2,:))/rdd2;
    rnn  = norm(Sys.Icoord(1,:)-Sys.Icoord(2,:));
    NN  = yI*yI*mu0*hbar/(4*pi*(rdd1*1e-10)^3)/(2*pi*1e6);
    uvnn = (Sys.Icoord(1,:)-Sys.Icoord(2,:))/rnn;
    cv   = ori.vecs';
    cs1  = sum(cv.*uv1);
    cs2  = sum(cv.*uv2);
    csnn = sum(cv.*uvnn);

    A1   = (3*cs1^2-1)*wdd1;
    A2   = (3*cs2^2-1)*wdd2;
    wnn  = (3*csnn^2-1)*NN;
    wZQ  = 0.5*sqrt((A1-A2)^2 + wnn^2);
    mod  = (((A1-A2)^2)*(wnn^2))/(16*(wZQ^4));
    Vori = 1 - 1.5*mod - 0.5*mod*cos(2*pi*wZQ*2*tau) + 2*mod*cos(2*pi*wZQ*tau);
end
