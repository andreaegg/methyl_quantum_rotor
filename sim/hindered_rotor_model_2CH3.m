function [nut,E] = hindered_rotor_model_2CH3(V3,W3,K,B)

% Calculate Energy Eigenvalues
E = zeros((2*K+1)^2,length(V3));
for k = 1:length(V3)
    currpot = V3(k);
    H = rotorhamiltonian_2CH3(currpot,W3,K,B);
    tic
    En = eig(H);
    En = sort(En);
    toc
    E(:,k) = En;
end

% Calculate Tunnel Frequency
for k = 1:length(V3)
    nut(k) = E(2,k) - E(1,k);
end


end