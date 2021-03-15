function K = tunnelkernel(vt,t,B,NO,CH3,knots)

ye = gfree*bmagn/hbar;
gyro = nucgval('1H')*nmagn/hbar;


if nargin < 3
    B = ones(size(t));
    NO = [];
    CH3 = [];
    knots = [];
end

if isempty(B)
    B = ones(size(t));
end

if nargin < 4
    NO = [];
    CH3 = [];
    knots = [];
end
if nargin < 6
    knots = 31;
end

if ~isempty(knots)
    [vecs,weights] = sphgrid('Ci',knots,'c');
    nori = length(weights); % number of orientations
end

model = @(a,omt,x)(2/3)*(12 + 2*cos(x*(a(1)-a(2))) - cos(x*(a(1)-a(2)-omt)) - 2*cos(x*(a(1)-a(3))).*(-1 + cos(x*omt)) ...
    - 2*cos(x*(a(2)-a(3))).*(-1 + cos(x*omt)) + 6*cos(x*omt) - cos(x*(a(1)-a(2)+omt)));

K = zeros(numel(vt),numel(t));

for m = 1:length(vt)
    currsig = 0;
    for ori = 1:nori
        cvec = vecs(:,ori)';
        distance = vecnorm((CH3 - NO),2,2);
        uvec = (CH3 - NO)./distance;
        wdd = ye*gyro*mu0*hbar./(4*pi*(distance*1e-10).^3)/1e6;
        ct = sum(cvec.*uvec,2);
        a0 = (3*ct.^2-1).*wdd;
        currsig = currsig + weights(ori)*model(a0,2*pi*vt(m),t); % coupling constant and tunnelfreq in MHz*rad
    end
    K(m,:) = currsig./max(currsig);
end

K = K.*B;

end




