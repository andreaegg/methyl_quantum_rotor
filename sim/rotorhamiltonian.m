function ham = rotorhamiltonian(V3,m,n,B)

if length(m) ~= length(n)
    error('Matrix must be symetric');
end

if nargin < 4
    B = 7.589357490029396; % default = 1H
end

ham = zeros(length(m),length(n));

for i = 1:length(m)
    for j = 1:length(n)
        Hij = (n(j)^2)*B*eq(m(i),n(j));
        Vij = (V3/4)*(2*eq(m(i),n(j))-eq(n(j)-m(i)+3,0)-eq(n(j)-m(i)-3,0));
        ham(i,j) = Hij + Vij;
    end
end

end