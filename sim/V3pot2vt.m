
function [vt,Potential,Energy] = V3pot2vt(V3,nreigval,rotor)

if nargin < 2
    nreigval = 50;
end
if nargin < 3
    rotor = 'CH3';
end

switch rotor
    case 'CH3'
       B = 7.6;        % K % Moreno et al,Phys.Rev.B,1999,59,9
    case 'CD3'
       B = 3.8;        % K % Moreno et al,Phys.Rev.B,1999,59,9
end

m = -nreigval:1:nreigval;
n = -nreigval:1:nreigval;

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

vt.kHz = nut*boltzm/planck/1e3;
vt.MHz = nut*boltzm/planck/1e6;
vt.K   = nut;

Energy.K     = E;
Energy.kJ    = E*boltzm/1e3;
Energy.kJmol = E*molgas/1e3;

Potential.K      = V3;
Potential.kJ     = V3*boltzm/1e3;
Potential.kJmol  = V3*molgas/1e3;
end