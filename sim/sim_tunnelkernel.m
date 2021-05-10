function [K,T,nu] = sim_tunnelkernel(vt,experiment,Sys,Exp,Opt,rotor)

if nargin < 6
    rotor = 'CH3';
end

Opt.force_complete = true;
T = linspace(2*Exp.tau,2*Exp.tau +2*(Exp.npoints-1)*Exp.dt,Exp.npoints);

K.td = zeros(length(vt),length(T));
K.fd = zeros(length(vt),Opt.zerofilling);

switch experiment
    case '2pESEEM'
        for k = 1:length(vt)
            curr_vt = vt(k);
            Sys.vt  = curr_vt;
            if rotor == 'CH3'
                [~,signal] = ex_2pESEEM_mqr_CH3_parallel(Sys,Exp,Opt);
            elseif rotor == 'CD3'
                [~,signal] = ex_2pESEEM_mqr_CD3_parallel(Sys,Exp,Opt);
            else
                error('Only CH3 and CD3 groups can be considered')
            end
            K.td(k,:) = signal;
        end
    case '3pESEEM'
        for k = 1:length(vt)
            curr_vt = vt(k);
            Sys.vt  = curr_vt;
            if rotor == 'CH3'
                [T,signal,nu,spc] = ex_3pESEEM_mqr_CH3_parallel(Sys,Exp,Opt);
            elseif rotor == 'CD3'
                [T,signal,nu,spc] = ex_3pESEEM_mqr_CD3_parallel(Sys,Exp,Opt);
            else
                error('Only CH3 and CD3 groups can be considered')
            end
            K.td(k,:) = signal;
            K.fd(k,:) = spc;
        end

end
end