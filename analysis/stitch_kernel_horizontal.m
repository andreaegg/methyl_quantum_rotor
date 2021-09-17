function [] = stitch_kernel_horizontal(input,parts)

output = input(1:end-1);

position = strfind(input,'CH3');
methyl   = input(position:position+3);

for k = 1:parts
    
    part     = strcat('part',num2str(k));
    filename = strcat(input,part);
    currdata = load(filename);
    switch methyl
        case 'CH3a'
            K = currdata.Ka;
        case 'CH3b'
            K = currdata.Kb;
        case 'CH3c'
            K = currdata.Kc;
        case 'CH3d'
            K = currdata.Kd;
    end
    idx      = find(K(:,1) ~= 0);
    start    = idx(1);
    kernel(start:start+length(idx)-1,:)   = K(start:start+length(idx)-1,:);
    
end



save(output,'kernel')
end