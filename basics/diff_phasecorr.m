function difference = diff_phasecorr(trace,phi)

difference = abs(imag(trace*exp(-1i*phi)));
difference = sum(difference);
end