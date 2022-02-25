function [trcorr,phi] = phasecorr(trace,idx)
if nargin < 2
    phi0 = atan2(imag(trace(end)),real(trace(end))); % start phase angle
    phi  = fminsearch(@(x)rms_phasecorr(trace,x),phi0);
    trcorr = trace*exp(-1i*phi);
    phi    = phi*180/pi;
    [~,idxmax] = max(abs(real(trcorr)));
    neg = trcorr(idxmax);
    if neg < 0
        trcorr = trcorr*exp(-1i*pi);
        phi = phi+180;
    end
else
    phi0 = atan2(imag(trace(idx(1))),real(trace(idx(1)))); % start phase angle
    currtrace = trace(idx);
    phi  = fminsearch(@(x)rms_phasecorr(currtrace,x),phi0);
    trcorr = trace*exp(-1i*phi);
    phi    = phi*180/pi;
    [~,idxmax] = max(abs(real(trcorr)));
    neg = trcorr(idxmax);
    if neg < 0
        trcorr = trcorr*exp(-1i*pi);
        phi = phi+180;
    end
end
end