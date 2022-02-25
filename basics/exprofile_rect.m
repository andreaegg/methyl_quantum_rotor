function [nu,spc] = exprofile_rect(tpulse)

t = linspace(-25*tpulse,25*tpulse,10001); % ns
t = t*1e-9;
box = zeros(size(t));
idxmin = find(t == -tpulse*1e-9/2);
idxmax = find(t == tpulse*1e-9/2);

box(idxmin:idxmax) = 1;

spc = fftshift(fft(box));
spc = spc/max(spc);

dt = mean(diff(t));
nu=0:1/(length(t)*dt):1/dt;
nu=nu(1:length(spc));
nu=fftshift(nu);
nu(nu>=(1/(2*dt)))=nu(nu>=(1/(2*dt)))-1/dt;
% nu = -1/(2*dt):1/(length(t)*dt):1/(2*dt);
% nu = nu(1:length(spc));
nu = nu/1e6; % MHz
end