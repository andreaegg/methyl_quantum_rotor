function [corrint,baseline] = baselinecorr(B,int,points,side,degree,mode)

if nargin < 3
    points = 100;
end
if nargin < 4
    side = 'right';
end
if nargin < 5
    degree = 1;
end
if nargin < 6
    mode = 'subtraction';
end

int = real(int);

switch side
    case 'left'
        param = polyfit(B(1:points),int(1:points),degree);
        baseline = polyval(param,B);
    case 'right'
        param = polyfit(B(end-points+1:end),int(end-points+1:end),degree);
        baseline = polyval(param,B);
    case 'both'
        param = polyfit([B(1:points) B(end-points+1:end)],[int(1:points) int(end-points+1:end)],degree);
        baseline = polyval(param,B);
end

switch mode
    case 'division'
        corrint = int./baseline;
    case 'subtraction'
        corrint = int - baseline;
end
end