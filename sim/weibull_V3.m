function pdf = weibull_V3(V,lambda,k)

if nargin < 2
    lambda = 1;
    k = 1;
end

if nargin < 3
    k = 1;
end

if (all(V > 0))
    P = (k/lambda).*((V/lambda).^(k-1)).*exp(-(V/lambda).^k);    
    pdf = P/sum(P);
else
    P = zeros(size(V));
    pdf = P;
end


end