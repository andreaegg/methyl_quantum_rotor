function hamtot = rotorhamiltonian_2CH3Khazaei(V3,W3,K,B)

if nargin < 4
    B = 7.589357490029396; % default = 1H
end

m = -K:1:K;
n = -K:1:K;
ham1 = zeros(length(m),length(n));
ham2 = zeros(length(m),length(n));
I   = eye(2*K+1);

% uncoupled rotors
for i = 1:length(m)
    for j = 1:length(n)
        C(i,j)   = (m(i)^2)*eq(m(i),n(j));
        Ip3(i,j) = 1*(eq(m(i)+3,n(j)));
        Im3(i,j) = 1*(eq(m(i)-3,n(j)));
    end
end

% full rotational Hamiltonian
ham1 = B*C + (V3/4)*(Ip3+Im3);
ham2 = B*C + (V3/4)*(Ip3+Im3);
hamtot = kron(ham1,I) + kron(I,ham2);
if abs(W3)>0
    Wij    = kron(Im3,Ip3) + kron(Ip3,Im3);
    hamtot = hamtot + (W3/2)*Wij;
end

end