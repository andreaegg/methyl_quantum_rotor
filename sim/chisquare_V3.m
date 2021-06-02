function pdf = chisquare_V3(V,k)

if nargin < 2
    k = 1;
end


if (all(V > 0))
    P = V.^(k/2).*exp(-(V/2))./(2.^(k/2).*gamma(k/2));    
    pdf = P/sum(P);
else
    P = zeros(size(V));
    pdf = P;
end


end