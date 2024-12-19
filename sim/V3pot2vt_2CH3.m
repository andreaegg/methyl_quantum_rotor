
function [vt,Potential,Energy] = V3pot2vt_2CH3(V3,W3,nreigval,rotor)

if nargin < 3
    nreigval = 22;
end
if nargin < 4
    rotor = 'CH3';
end

switch rotor
    case 'CH3'
       B = 7.6;        % K % Moreno et al,Phys.Rev.B,1999,59,9
    case 'CD3'
       B = 3.8;        % K % Moreno et al,Phys.Rev.B,1999,59,9
    case 'NH3'
       B = 8.9181;
end

K = nreigval;

[nut,E] = hindered_rotor_model_2CH3(V3,W3,K,B);

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