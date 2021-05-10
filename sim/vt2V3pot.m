function [V3,energy,vt] = vt2V3pot(vt,nreigval,rotor,V3start)

% Input:
% vt : in MHz

options = optimset('MaxFunEvals',1000,'MaxIter',1000,'Display','off');

if nargin < 2
    nreigval = 50;
end
if nargin < 3
    rotor = 'CH3';
end
if nargin < 4
    V3start = 0;
end

switch rotor
    case 'CH3'
        B = 7.6;        % K % Moreno et al,Phys.Rev.B,1999,59,9
    case 'CD3'
        B = 3.8;        % K % Moreno et al,Phys.Rev.B,1999,59,9
end

m = -nreigval:1:nreigval;
n = -nreigval:1:nreigval;
vt = vt*1e6; % comparison in Hz

for k = 1:length(vt)
    currvt = vt(k);
    V3opt = fminsearch(@(x)diff_vtfit(x,currvt,m,n,B),V3start,options);
    currV3(k) = V3opt;
end



V3.K     = currV3;
V3.kJ    = currV3*boltzm/1e3;
V3.kJmol = currV3*molgas/1e3;

[vt,energy] = hindered_rotor_model(currV3,m,n,B);


end