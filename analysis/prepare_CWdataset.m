function data = prepare_CWdataset(input1,input2,input3,norm)

if nargin < 2
    if isstruct(input1)
        % adjust B-axis
        g     = (planck*input1.param.MWFQ)./(bmagn*input1.B*1e-4);
        stfrq = 9.372*1e9;
        data.B = planck*stfrq*1e4./(g*bmagn);
        data.param = input1.param;
    end
elseif nargin < 5
    % adjust B-axis
    g     = (planck*input3.MWFQ)./(bmagn*input1*1e-4);
    stfrq = 9.372511*1e9;
    data.B = planck*stfrq*1e4./(g*bmagn);
    data.param = input3;
end
if norm
    data.spc = input2/max(input2);
end
end