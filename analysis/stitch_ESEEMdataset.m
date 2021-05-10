function [t,data,param] = stitch_ESEEMdataset(file1,file2,phase,offset,baseline,points,plotflag)

if nargin < 3
    phase = false;
end
if nargin < 4
    offset = false;
end
if nargin < 5
    baseline = false;
end
if nargin < 6
    points = 100;
end
if nargin < 7
    plotflag = false;
end

% load data and validation
if iscell(file1)
    t1    = file1{1};
    t1    = t1(:)/1000;
    data1 = file1{2};
    data1 = data1(:);
    t2    = file2{1};
    t2    = t2(:)/1000;
    data2 = file2{2};
    data2 = data2(:);
end
start1 = data1;
start2 = data2;

% optional phase correction
if phase
    [data1,phi1] = phasecorr(data1);
    [data2,phi2] = phasecorr(data2);
end

% find overlap
idx = find(t2 == t1(end));

% optional offset correction real signal
if offset
    diff = real(data1(end)) - real(data2(idx));
    data2  = data2 + diff;
end

% stitching
t    = [t1;t2(idx+1:end)];
data = [data1;data2(idx+1:end)];

% optional baseline correction real signal
if baseline
    [data,base] = baselinecorr(t,data,points,'right',1,'subtraction');
end

if plotflag
    figure(1),clf
    plot(t1,real(start1),'k')
    hold on;
    plot(t1,imag(start1),'k','HandleVisibility','off')
    plot(t2,real(start2),'k','HandleVisibility','off')
    plot(t2,imag(start2),'k','HandleVisibility','off')
    plot(t,real(data),'b')
    plot(t,imag(data),'b','HandleVisibility','off')
    xlim([0 max(t)])
    xlabel('T [us]')
    ylabel('intensity [a.u.]')
    legend({'raw','stitched,corrected'})
end

param.phasecorr    = phase;
if phase
    param.phase    = [phi1,phi2];
end
param.offsetcorr   = offset;
if offset
    param.offset   = diff;
end
param.baselinecorr = baseline;
if baseline
    param.baseline = base;
end

end