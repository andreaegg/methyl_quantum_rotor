
function xyz = dihedral_rotation(xyz0,dihedral,ethyl,angle,move)
global vec
% correct origin to first bond atom
origin = xyz0(dihedral(1),:);
xyz    = xyz0 - origin;

% rotation of molecule that dihedral bond corresponds to x-axis
vec = xyz(dihedral(2),:);
a0 = [0,0,0];
options = optimoptions('fsolve','Algorithm','levenberg-marquardt');
a = fsolve(@rotation_bond2newx,a0);
R = rotation_matrix(a(1),a(2),a(3));
xyz = (R*xyz')';

% rotation around dihedral angle of bond
[~,Rx] = rotation_matrix(0,0,angle);
xyz_ethyl = Rx*(xyz(ethyl,:)');
xyz_move  = Rx*(xyz(move,:)');

xyz(ethyl,:) = xyz_ethyl';
xyz(move,:)  = xyz_move';


function F = rotation_bond2newx(a)
    % alpha a(1) = Rz
    % beta  a(2) = Ry
    % gamma a(3) = Rx
    F(1) = cos(a(1))*cos(a(2))*vec(1) + (cos(a(1))*sin(a(2))*sin(a(3))-sin(a(1))*cos(a(3)))*vec(2) + (cos(a(1))*sin(a(2))*cos(a(3))+sin(a(1))*sin(a(3)))*vec(3)-norm(vec);
    F(2) = cos(a(2))*sin(a(3))*vec(1) + (sin(a(1))*sin(a(2))*sin(a(3))+cos(a(1))*cos(a(3)))*vec(2) + (sin(a(1))*sin(a(2))*cos(a(3))-cos(a(1))*sin(a(3)))*vec(3);
    F(3) = -sin(a(2))*vec(1) + sin(a(1))*cos(a(2))*vec(2) + cos(a(1))*cos(a(2))*vec(3);
end

function [Rtot,Rx,Ry,Rz] = rotation_matrix(alpha,beta,gamma)
    Rx = [1 0 0;0 cos(gamma) -sin(gamma);0 sin(gamma) cos(gamma)];
    Ry = [cos(beta) 0 sin(beta);0 1 0;-sin(beta) 0 cos(beta)];
    Rz = [cos(alpha) -sin(alpha) 0;sin(alpha) cos(alpha) 0;0 0 1;];
    Rtot = Rz*Ry*Rx;
end


end