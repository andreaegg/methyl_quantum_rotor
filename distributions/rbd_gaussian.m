function P = rbd_gaussian(V,V0,sigma)

x = (V-V0)/(sigma*sqrt(2));
P = (1/(sigma*sqrt(2*pi)))*exp(-(x.^2));
P = P/trapz(P)/mean(diff(V));
end

