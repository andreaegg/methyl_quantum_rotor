function rms = rms_phasecorr(trace,phi)
imagtr = imag(trace*exp(-1i*phi));
rms    = sqrt((1/numel(trace))*sum(imagtr.*imagtr));
end