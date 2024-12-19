function hamtot = rotorhamiltonian_2CH3(V3,W3,K,B)

if nargin < 4
    B = 7.589357490029396; % default = 1H
end

m = -K:1:K;
n = -K:1:K;
ham = zeros(length(m),length(n));
I   = eye(2*K+1);

% uncoupled rotors
for i = 1:length(m)
    for j = 1:length(n)
        Hij = (n(j)^2)*B*eq(m(i),n(j));
        Vij = (V3/4)*(2*eq(m(i),n(j)) + eq(m(i)+3,n(j)) + eq(m(i),n(j)+3));
        ham(i,j) = Hij + Vij;
    end
end

% coupled rotors
if abs(W3)>0
    Wij_1a = zeros(length(m),length(n));
    Wij_2a = zeros(length(m),length(n));
    Wij_1b = zeros(length(m),length(n));
    Wij_2b = zeros(length(m),length(n));

    for i = 1:length(m)
        for j = 1:length(n)
        Wij_1a(i,j) = eq(m(i)-3,n(j));
        Wij_2a(i,j) = eq(m(i)+3,n(j));

        Wij_1b(i,j) = eq(m(i)+3,n(j));
        Wij_2b(i,j) = eq(m(i)-3,n(j));
        end
    end
    Wij = kron(Wij_1a,Wij_2a) + kron(Wij_1b,Wij_2b);
end

% full rotational Hamiltonian
hamtot = kron(ham,I) + kron(I,ham);
if abs(W3)>0
    hamtot = hamtot + (W3/2)*Wij;
end

end