
function [] = stitch_kernel_horizontal(input,parts)

output = input(1:end-1);

position = strfind(input,'CH3');
methyl   = input(position:position+3);

for k = 1:parts
    
    part     = strcat('part',num2str(k));
    filename = strcat(input,part);
    currdata = load(filename);
    V3 = currdata.V3;
    T  = currdata.T;
    vt = currdata.vt;
    Exp = currdata.Exp;
    Sys = currdata.Sys;
    Opt = currdata.Opt;
    weights = currdata.weights;
    
    switch methyl
        case 'CH3a'
            K = currdata.Ka;
            idx      = find(K(:,1) ~= 0);
            start    = idx(1);
            Ka(start:start+length(idx)-1,:,:)   = K(start:start+length(idx)-1,:,:);
            save(output,'-v7.3','Ka','weights','V3','T','vt','Exp','Sys','Opt')
        case 'CH3b'
            K = currdata.Kb;
            idx      = find(K(:,1) ~= 0);
            start    = idx(1);
            Kb(start:start+length(idx)-1,:,:)   = K(start:start+length(idx)-1,:,:);
            save(output,'-v7.3','Kb','weights','V3','T','vt','Exp','Sys','Opt')
        case 'CH3c'
            K = currdata.Kc;
            idx      = find(K(:,1) ~= 0);
            start    = idx(1);
            Kc(start:start+length(idx)-1,:,:)   = K(start:start+length(idx)-1,:,:);
            save(output,'-v7.3','Kc','weights','V3','T','vt','Exp','Sys','Opt')
        case 'CH3d'
            K = currdata.Kd;
            idx      = find(K(:,1) ~= 0);
            start    = idx(1);
            Kd(start:start+length(idx)-1,:,:)   = K(start:start+length(idx)-1,:,:);
            save(output,'-v7.3','Kd','weights','V3','T','vt','Exp','Sys','Opt')
    end
    
end
end