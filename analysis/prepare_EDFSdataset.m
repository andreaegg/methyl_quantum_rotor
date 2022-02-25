function data = prepare_EDFSdataset(input1,input2,input3,norm)

if nargin < 2
    if isstruct(input1)
        % phase correction
        data.spc = phasecorr(input1.spc);
        % adjust B-axis
        g     = (planck*input1.param.MWFQ)./(bmagn*input1.B*1e-4);
        stfrq = 34.5*1e9;
        data.B = planck*stfrq*1e3./(g*bmagn);
        data.param = input1.param;
    end
elseif nargin < 4
    % phase correction
    data.spc = phasecorr(input2);
    % adjust B-axis
    g     = (planck*input3.MWFQ)./(bmagn*input1*1e-4);
    stfrq = 34.5*1e9;
    data.B = planck*stfrq*1e3./(g*bmagn);
    data.param = input3;
end
if norm
    data.spc = data.spc/max(data.spc);
end
end