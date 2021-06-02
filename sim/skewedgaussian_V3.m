function pdf = skewedgaussian_V3(V,V0,stdV,alpha,factor)

if nargin < 2
    V0 = mean(V);
    stdV = std(V)/2;
    alpha = 0;
    factor = 1;
end

if nargin < 3
    stdV = std(V)/2;
    alpha = 0;
    factor = 1;
end

if nargin < 4
    alpha = 0;
    factor = 1;
end

if nargin < 5
    factor = 1;
end

P = zeros(size(V));

gaussian     = @(V,V0,stdV) (1/(stdV*sqrt(2*pi)))*exp(-(V-V0).^2/(2*stdV^2));

if size(V0) == size(stdV)
    if size(V0) == size(factor)
        if size(stdV) == size(factor)
            for k = 1:length(V0)
                skew      = normcdf(alpha*(V-V0)/stdV);
                P = P + 2/stdV*gaussian(V,V0,stdV).*skew;
            end
        end
    end
end

pdf = P/sum(P);

end