function [V,freq] = ESEEM2ptheoretical_mckracken(t,nuc,r,ori,points,field,Sys)
for n = 1:length(r)
    yI  = nucgval(nuc)*nmagn/hbar;
    ye  = Sys.g*bmagn/hbar;
    wdd = ye*yI*mu0*hbar/(4*pi*(r(n)*1e-10)^3)/(2*pi*1e6);
    Vr = zeros(1,points);
    for m = 1:length(ori.weights)
        cv  = ori.vecs(:,m)';
        cs  = sum(cv.*ori.uv(n,:));
        A   = (3*cs^2-1)*wdd;
        B   = 3*cs*sqrt(1-cs^2)*wdd;
        wI  = -yI*field/2/pi/1e6;
        wa(n,m) = sqrt((wI-A/2)^2+(B^2)/4);
        wb(n,m) = sqrt((wI+A/2)^2+(B^2)/4);
        wp(n,m)   = wa(n,m) + wb(n,m);
        wm(n,m)   = wa(n,m) - wb(n,m);
        mod = (wI*B/(wa(n,m)*wb(n,m)))^2;
        Vcurr = 1-(mod/2)*(1-cos(2*pi*wa(n,m)*t)-cos(2*pi*wb(n,m)*t)+0.5*cos(2*pi*wm(n,m)*t)+0.5*cos(2*pi*wp(n,m)*t));
        Vr = Vr + ori.weights(m)*Vcurr;
    end
    freq = {wa,wb};
    V = Vr/sum(ori.weights);
end
end