function [K,T,nu] = sim_tunnelkernel(vt,experiment,Sys,Exp,Opt,rotor)

if nargin < 6
    rotor = 'CH3';
end

% intialization
Opt.force_complete = true;
nu = [];

K.td = zeros(length(vt),Exp.npoints);

% simulation
switch experiment
    case '2pESEEM'
        for k = 1:length(vt)
            curr_vt = vt(k);
            Sys.vt  = curr_vt;
            if strcmp(rotor,'CH3')
                [T,signal] = ex_2pESEEM_mqr_CH3_parallel(Sys,Exp,Opt);
            elseif strcmp(rotor,'CD3')
                [T,signal] = ex_2pESEEM_mqr_CD3_parallel(Sys,Exp,Opt);
            else
                error('Only CH3 and CD3 groups can be considered')
            end
            K.td(k,:) = signal;
        end
        
    case '3pESEEM'
        for k = 1:length(vt)
            curr_vt = vt(k);
            Sys.vt  = curr_vt;
            if strcmp(rotor,'CH3')
                [T,signal] = ex_3pESEEM_mqr_CH3_parallel(Sys,Exp,Opt);
            elseif strcmp(rotor,'CD3')
                [T,signal] = ex_3pESEEM_mqr_CD3_parallel(Sys,Exp,Opt);
            else
                error('Only CH3 and CD3 groups can be considered')
            end
            K.td(k,:) = signal;
        end
        
    case 'DEFENCE'
        for k = 1:length(vt)
            curr_vt = vt(k);
            Sys.vt  = curr_vt;
            if strcmp(rotor,'CH3')
                [T,signal] = ex_DEFENCE_mqr_CH3_parallel(Sys,Exp,Opt);
            else
                error('Only CH3 groups can be considered')
            end
            K.td(k,:) = signal;
        end
        
    case '5pESEEM'
        for k = 1:length(vt)
            curr_vt = vt(k);
            Sys.vt  = curr_vt;
            if strcmp(rotor,'CH3')
                [T,signal] = ex_5pESEEM_mqr_CH3_parallel(Sys,Exp,Opt);
            else
                error('Only CH3 groups can be considered')
            end
            K.td(k,:) = signal;
        end
        
    case 'CPDD2'
        for k = 1:length(vt)
            curr_vt = vt(k);
            Sys.vt  = curr_vt;
            if strcmp(rotor,'CH3')
                [T,signal] = ex_CPDD2_mqr_CH3_parallel(Sys,Exp,Opt);
            else
                error('Only CH3 groups can be considered')
            end
            K.td(k,:) = signal;
        end
    case 'CPDD3'
        for k = 1:length(vt)
            curr_vt = vt(k);
            Sys.vt  = curr_vt;
            if strcmp(rotor,'CH3')
                [T,signal] = ex_CPDD3_mqr_CH3_parallel(Sys,Exp,Opt);
            else
                error('Only CH3 groups can be considered')
            end
            K.td(k,:) = signal;
        end
    case 'CPDD4'
        for k = 1:length(vt)
            curr_vt = vt(k);
            Sys.vt  = curr_vt;
            if strcmp(rotor,'CH3')
                [T,signal] = ex_CPDD4_mqr_CH3_parallel(Sys,Exp,Opt);
            else
                error('Only CH3 groups can be considered')
            end
            K.td(k,:) = signal;
        end
end
end