function pdf = gaussian_V3(V,V0,stdV,factor)

if nargin < 2
    V0 = mean(V);
    stdV = std(V)/2;
    factor = 1;
end

if nargin < 3
    stdV = std(V)/2;
    factor = 1;
end

if nargin < 4
    factor = 1;
end

P = zeros(size(V));

if size(V0) == size(stdV)
    if size(V0) == size(factor)
        if size(stdV) == size(factor)
            for k = 1:length(V0)
                P = P + (factor(k)/(stdV(k)*sqrt(2*pi)))*exp(-(V-V0(k)).^2/(2*stdV(k)^2));
            end
        end
    end
end

pdf = P/sum(P);

end