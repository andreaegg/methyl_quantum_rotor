function [propagator,propagatorr] = gen_propagators_real(tpulse,pulses,ham)

for k = 1:length(pulses)
    propagator{k}  = expm(-1i*2*pi*tpulse*(pulses{k}+ham));
    propagatorr{k} = expm( 1i*2*pi*tpulse*(pulses{k}+ham));
end

end