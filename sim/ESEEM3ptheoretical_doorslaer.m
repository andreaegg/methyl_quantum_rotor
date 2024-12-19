function [V,freq] = ESEEM3ptheoretical_doorslaer(t,tau,nuc,r,ori,points,field)

Vr = zeros(1,points);
for m = 1:length(ori.weights)
    Voria = ones(1,points);
    Vorib = ones(1,points);
    for n = 1:length(r)
        yI  = nucgval(nuc)*nmagn/hbar;
        ye  = gfree*bmagn/hbar;
        wdd = ye*yI*mu0*hbar/(4*pi*(r(n)*1e-10)^3)/(2*pi*1e6);
        cv  = ori.vecs(:,m)';
        cs  = sum(cv.*ori.uv(n,:));
        A   = (3*cs^2-1)*wdd;
        B   = 3*cs*sqrt(1-cs^2)*wdd;
        wI  = -yI*field/2/pi/1e6;
        na  = atan(-B/(2*wI+A));
        nb  = atan(-B/(A-2*wI));
        eta = (na - nb)/2;
        wa(n,m)  = abs((wI + A/2)*cos(na)-(B/2)*sin(na));
        wb(n,m)  = abs((wI - A/2)*cos(nb)+(B/2)*sin(nb));
        mod   = sin(2*eta)^2;
        Vcurra = (1/2)*(1 - (mod/2)*(1 - cos(2*pi*wb(n,m)*tau)).*(1 - cos(2*pi*wa(n,m)*(t+tau))));
        Vcurrb = (1/2)*(1 - (mod/2)*(1 - cos(2*pi*wa(n,m)*tau)).*(1 - cos(2*pi*wb(n,m)*(t+tau))));
        Voria  = Voria.*Vcurra;
        Vorib  = Vorib.*Vcurrb;
    end
    freq = {wa,wb};
    Vori = Voria + Vorib;
    Vr = Vr + ori.weights(m)*Vori;
end
V = Vr/sum(ori.weights);
end