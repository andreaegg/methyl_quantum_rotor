function [K,T] = sim_tunnelkernel(vt,experiment,Sys,Exp,Opt,rotor)

if nargin < 6
    rotor = 'CH3';
end

if ismember(experiment,["3pSIFTER","4pSIFTER"])
    K = [];
else
    % intialization
    K = zeros(length(vt),Exp.npoints);
end

% simulation
switch experiment
    case "2pESEEM"
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
            K(k,:) = phasecorr(signal);
        end
        
    case "3pESEEM"
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
            K(k,:) = phasecorr(signal);
        end
        
    case "DEFENCE"
        for k = 1:length(vt)
            curr_vt = vt(k);
            Sys.vt  = curr_vt;
            if strcmp(rotor,'CH3')
                [T,signal] = ex_DEFENCE_mqr_CH3_parallel(Sys,Exp,Opt);
            else
                error('Only CH3 groups can be considered')
            end
            K(k,:) = phasecorr(signal);
        end
        
    case "5pESEEM"
        for k = 1:length(vt)
            curr_vt = vt(k);
            Sys.vt  = curr_vt;
            if strcmp(rotor,'CH3')
                [T,signal] = ex_5pESEEM_mqr_CH3_parallel(Sys,Exp,Opt);
            else
                error('Only CH3 groups can be considered')
            end
            K(k,:) = phasecorr(signal);
        end
    case "3pSIFTER"
        for k = 1:length(vt)
            curr_vt = vt(k);
            Sys.vt  = curr_vt;
            if strcmp(rotor,'CH3')
                [T,signal] = ex_3pSIFTER_mqr_CH3_parallel(Sys,Exp,Opt);
            else
                error('Only CH3 groups can be considered')
            end
            K(k,:) = phasecorr(signal);
        end
    case "4pSIFTER"
        for k = 1:length(vt)
            curr_vt = vt(k);
            Sys.vt  = curr_vt;
            if strcmp(rotor,'CH3')
                [T,signal] = ex_4pSIFTER_mqr_CH3_parallel(Sys,Exp,Opt);
            else
                error('Only CH3 groups can be considered')
            end
            K(k,:) = phasecorr(signal);
        end
end
end