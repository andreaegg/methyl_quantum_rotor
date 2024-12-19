function ham = rotorhamiltonian_Khazaei(V3,K,B)

if nargin < 4
    B = 7.589357490029396; % default = 1H
end

m = -K:1:K;
n = -K:1:K;

% uncoupled rotors
for i = 1:length(m)
    for j = 1:length(n)
        C(i,j)   = (m(i)^2)*eq(m(i),n(j));
        Ip3(i,j) = 1*(eq(m(i)+3,n(j)));
        Im3(i,j) = 1*(eq(m(i)-3,n(j)));
    end
end

% full rotational Hamiltonian
ham = B*C + (V3/4)*(Ip3+Im3);
end