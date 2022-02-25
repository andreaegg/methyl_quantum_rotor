function [rgb] = atomcolors(atom)

if nargin < 1
    atom = "H";
end

switch atom
    case "H"
        rgb = colorss('cadet grey');
    case "C"
        rgb = colorss('black');
    case "N"
        rgb = colorss('blue');
    case "O"
        rgb = colorss('red');
    case "Cu"
        rgb = colorss('brown (traditional)');       
end

end