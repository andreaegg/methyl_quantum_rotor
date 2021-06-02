function pdf = powerlognormal_V3(V,V0,stdV,p)

if nargin < 2
    V0 = mean(V);
    stdV = std(V)/2;
    p = 1;
end

if nargin < 3
    stdV = std(V)/2;
    p = 1;
end

gaussian = normpdf((log(V)-V0)/stdV);
skew     = (normcdf(-(log(V)-V0)/stdV)).^(p-1);

P = (p./(V*stdV)).*gaussian.*skew;

pdf = P/sum(P);

end