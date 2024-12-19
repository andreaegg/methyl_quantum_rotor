function Ea = V3pot2Ea(V3,nreigval,rotor)

if nargin < 2
    nreigval = 51;
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

[~,E] = hindered_rotor_model(V3,m,n,B);

Ea.r0 = V3 - min(E);
Ea.r1 = V3 - E(4);

end