function means = mean_phasecorr(trace,phi)

means = abs(imag(trace*exp(-1i*phi)));
means = mean(means);
end