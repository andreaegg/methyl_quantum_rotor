function [nut,E] = hindered_rotor_model(V3,m,n,B)


% Calculate Energy Eigenvalues
E = zeros(length(m),length(V3));
for k = 1:length(V3)
    currpot = V3(k);
    H = rotorhamiltonian(currpot,m,n,B);
    En = eig(H);
    En = sort(En);
    E(:,k) = En;
end

% Calculate Tunnel Frequency
for k = 1:length(V3)
    nut(k) = E(2,k) - E(1,k);
end


end