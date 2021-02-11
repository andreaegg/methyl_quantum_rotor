function dist = distribution_V3(V,V0,stdV)

if nargin < 2
    V0 = mean(V);
    stdV = std(V)/2;
end

P = (1/(stdV*sqrt(2*pi)))*exp(-(V-V0).^2/(2*stdV^2));
dist.x   = V;
dist.pdf = P/sum(P);

end