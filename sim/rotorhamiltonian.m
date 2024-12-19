function ham = rotorhamiltonian(V3,K,B)

if nargin < 3
    B = 7.589357490029396; % default = 1H
end

m = -K:1:K;
n = -K:1:K;
ham = zeros(length(m),length(n));

for i = 1:length(m)
    for j = 1:length(n)
        % Khazaei
        % Hij = B*(n(j)^2)*eq(m(i),n(j));
        % Vij = (V3/4)*(eq(m(i)+3,n(j)) + eq(m(i)-3,n(j)));
        % Andrea
        Hij = (n(j)^2)*B*eq(m(i),n(j));
        Vij = (V3/4)*(2*eq(m(i),n(j)) + eq(m(i)+3,n(j)) + eq(m(i),n(j)+3));
        ham(i,j) = Hij + Vij;
    end
end

end