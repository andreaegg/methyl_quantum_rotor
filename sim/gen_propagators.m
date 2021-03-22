function [propagator,propagatorr] = gen_propagators(tpulse,pulses)

for k = 1:length(pulses)
    propagator{k} = expm(-1i*2*pi*(pulses{k})*tpulse);
    propagatorr{k} = expm(1i*2*pi*(pulses{k})*tpulse);
end

end