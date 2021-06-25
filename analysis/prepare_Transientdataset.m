function data = prepare_Transientdataset(input)

% phase correction
data.signal = phasecorr(input.signal);

% normalization
data.signal = data.signal/max(data.signal);
data.t      = input.t;

end