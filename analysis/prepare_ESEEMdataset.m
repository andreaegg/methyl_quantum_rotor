function [t,data,param] = prepare_ESEEMdataset(file,phase,baseline,points,plotflag)

if nargin < 2
    phase = false;
end
if nargin < 3
    baseline = false;
end
if nargin < 4
    points = 100;
end
if nargin < 5
    plotflag = false;
end

t = file{1};
t = t(:)/1000;
data = file{2};
data = data(:);

start = data;

if phase
    [data,phi] = phasecorr(data);
end

if baseline
    [data,base] = baselinecorr(t,data,points,'right',1,'subtraction');
end

if plotflag
    figure(1),clf
    plot(t,real(start),'k')
    hold on;
    plot(t,imag(start),'k','HandleVisibility','off')
    plot(t,real(data),'b')
    plot(t,imag(data),'b','HandleVisibility','off')
    xlim([0 max(t)])
    xlabel('T [us]')
    ylabel('intensity [a.u.]')
    legend({'raw','corrected'})
end

param.phasecorr    = phase;
if phase
    param.phase    = phi;
end
param.baselinecorr = baseline;
if baseline
    param.baseline = base;
end

end