function [trcorr,phi] = phasecorr(trace,idx)
if nargin < 2
    phi0 = atan(imag(trace(end))/real(trace(end))); % start phase angle
    phi  = fminsearch(@(x)rms_phasecorr(trace,x),phi0);
    trcorr = trace*exp(-1i*phi);
else
    phi0 = atan(imag(trace(idx(1)))/real(trace(idx(1)))); % start phase angle
    currtrace = trace(idx);
    phi  = fminsearch(@(x)rms_phasecorr(currtrace,x),phi0);
    trcorr = trace*exp(-1i*phi);
end
end