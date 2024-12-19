function [sig_pc] = propagation_pc(sigma,propagator,propagatorr,addphase)

if isempty(addphase)
    sig_pc = propagator{1}*sigma*propagatorr{1};
else
    sig_pc = zeros(size(propagator{1}));
    for k = 1:length(propagator)
        switch addphase(k)
            case "+"
                sig_tmp = propagator{k}*sigma*propagatorr{k};
            case "-"
                sig_tmp = -(propagator{k}*sigma*propagatorr{k});
        end
        sig_pc = sig_pc + sig_tmp/length(propagator);
    end
end
end