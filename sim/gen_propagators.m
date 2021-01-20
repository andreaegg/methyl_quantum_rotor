function propagator = gen_propagators(currham,tpulse,pulses)

for k = 1:length(pulses)
    propagator{k} = expm(-1i*2*pi*(currham+pulses{k})*tpulse);
end

end