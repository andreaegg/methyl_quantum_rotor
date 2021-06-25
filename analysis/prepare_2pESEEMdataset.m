function data = prepare_2pESEEMdataset(input)

% phase correction
data.signal = phasecorr(input.signal);

% normalization
data.signal = data.signal/max(data.signal);
data.t      = input.t/1000;

end