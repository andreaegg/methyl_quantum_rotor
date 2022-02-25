function output = prepare_pc_presat2pESEEMdataset(input_px,input_mx)

% Define time-axis to us
if iscell(input_px.t)
    output.t{1} = input_px.t{1};
    output.t{2} = input_px.t{2};
else
    output.t = input_px.t;
end

% Add +x presat and -x presat traces
data_px = input_px.signal;
data_mx = input_mx.signal;
currsig = (data_px + data_mx)/2; % gets rid of moving echos
for n = 1:size(currsig,2)
    curr = phasecorr(currsig(:,n));
    output.signal(:,n) = curr;
end

end


