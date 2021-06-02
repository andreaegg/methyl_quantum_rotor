function pdf = grumbleminimum_V3(V,V0,beta)

if nargin < 2
    V0 = mean(V);
    beta = 1;
end

if nargin < 3
    beta = 1;
end

P = (1/beta)*exp((V-V0)/beta).*exp(-exp((V-V0)/beta));

pdf = P/sum(P);

end