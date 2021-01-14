function [sig_pc] = propagation_pc(sigma,propagator,addphase)

if isempty(addphase)
    addphase = ["+","-"];
end

sig_pc = zeros(size(propagator{1}));

if length(propagator) == 1
    sig_pc = propagator{1}*sigma*propagator{1}';
else
    for k = 1:length(propagator)
        switch addphase(k)
            case "-"
                sig_tmp = -propagator{k}*sigma*propagator{k}';
            case "+"
                sig_tmp = propagator{k}*sigma*propagator{k}';
        end
        sig_pc = sig_pc + sig_tmp/length(propagator);
    end
end

end