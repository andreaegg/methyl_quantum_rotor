function pdf = gammadistribution_V3(V,alpha,beta)

if nargin < 2
    alpha = 1;
    beta = 1;
end

if nargin < 3
    beta = 1;
end

if (all(V > 0))
    P = (beta.^alpha).*(V.^(alpha-1)).*exp(-beta*V)/gamma(alpha);    
    pdf = P/sum(P);
else
    P = zeros(size(V));
    pdf = P;
end


end