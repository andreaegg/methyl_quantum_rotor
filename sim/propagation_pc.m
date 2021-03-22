function [sig_pc] = propagation_pc(sigma,propagator,propagatorr,addphase)

if isempty(addphase)
    addphase = ["+","-"];
end

if ~iscell(propagator)
    sig_pc = propagator*sigma*propagatorr;
else
    sig_pc = zeros(size(propagator{1}));
    for k = 1:length(propagator)
        switch addphase(k)
            case "-"
                sig_tmp = -propagator{k}*sigma*propagatorr{k};
            case "+"
                sig_tmp = propagator{k}*sigma*propagatorr{k};
        end
        sig_pc = sig_pc + sig_tmp/length(propagator);
    end
end

end