function [connectivity,bondorder] = atomsconnectivity(atoms)

molecule = numel(atoms);

if (sum(atoms == "Cu") > 0)
    CuComplex = true;
    Nitroxide = false;
elseif (sum(atoms == "N") > 0)
    CuComplex = false;
    Nitroxide = true;
end

if CuComplex
    switch molecule
        case 36
            connectivity = [1 2; 1 7; 1 13;1 19; 1 27; 2 3;2 9;2 15;3 4;3 5;3 6;6 7;6 8;...
                            9 10;9 11;9 12;12 13;12 14;15 16;15 17;15 18;19 20;19 25;20 21;21 22;...
                            21 23;23 24;23 25;25 26;27 28;27 33;28 29;29 30;29 31;31 32;31 33;33 34;...
                            28 35;20 36];
            bondorder = [1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 2 1 1 1 2 1 1 1 ...
                1 2 1 1 1 2 1 1 1 1 2 1 1 1]; % not sure if true
        case 45
            connectivity = [1 2;1 3;1 7;1 18;2 11;2 12;2 13;3 4; 3 5;3 6;4 14;4 15;7 8;7 9;7 10;...
                           11 16;11 17;15 18;17 18;18 19;18 20;18 21;21 22;21 23;22 24;22 37;...
                           23 25;23 42;24 25;24 38;25 41;20 26;20 29;26 27;26 35;27 28;27 36;...
                           28 29;28 40;29 44;19 30;19 33;30 31;30 34;31 32;31 45;32 33;32 39;33 43];
            bondorder = [1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 2 1 1 1 1 1 1 1 2 2 ...
                1 1 1 1 1 1 1 1 2 2 1 1 1 1 1 1 2 2 1 1 1 1 1 1]; % not sure if true
    end
elseif Nitroxide
    switch molecule
%         case  % IAT
%             connectivity = [];
%             bondorder = [];
%         case  % IAP
%             connectivity = [];
%             bondorder = [];
%         case  % MALT
%             connectivity = [];
%             bondorder = [];
%         case  % MALP
%             connectivity = [];
%             bondorder = [];
%         case  % MTSL
%             connectivity = [];
%             bondorder = [];
        case 29 % TEMPO
            connectivity = [1 2;1 3;1 14;1 18;2 6;2 7;2 8;3 4;3 13;4 5;4 22;4 26;5 6;5 11;5 12;6 9;6 10;...
                            14 15;14 16;14 17;18 19;18 20;18 21;22 23;22 24;22 25;26 27;26 28;26 29];
            bondorder = [1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
    end
end
    
    
end