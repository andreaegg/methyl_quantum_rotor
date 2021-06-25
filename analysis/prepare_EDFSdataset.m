function data = prepare_EDFSdataset(input)

% phase correction
data.signal = phasecorr(input.signal);

% normalization
data.signal = data.signal/max(data.signal);

% adjust B-axis
g     = (planck*input.param.MWFQ)./(bmagn*input.t*1e-4);
stfrq = 34.5*1e9;
data.t = planck*stfrq*1e3./(g*bmagn);
data.param = input.param;
end