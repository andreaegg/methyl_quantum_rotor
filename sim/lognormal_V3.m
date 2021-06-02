function pdf = lognormal_V3(V,V0,stdV)

if nargin < 2
    V0 = mean(V);
    stdV = std(V)/2;
end

if nargin < 3
    stdV = std(V)/2;
end

P = (1./(V*stdV*sqrt(2*pi))).*exp(-(log(V)-V0).^2/(2*stdV^2));

pdf = P/sum(P);

end