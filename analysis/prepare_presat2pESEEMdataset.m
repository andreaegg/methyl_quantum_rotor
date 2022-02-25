function data = prepare_presat2pESEEMdataset(basename,variable,pc_flag,trec)
warning('off')

if nargin < 3
    pc_flag = zeros(length(variable),1);
end

told = 1;
for k = 1:length(variable)
    if variable(k) ~= ""
        if pc_flag(k)
            filename_p = strcat(basename,variable(k),"_+x");
            [plus.t,plus.signal] = eprload(convertStringsToChars(filename_p));
            filename_m = strcat(basename,variable(k),"_-x");
            [minus.t,minus.signal] = eprload(convertStringsToChars(filename_m));
            output = prepare_pc_presat2pESEEMdataset(plus,minus);
        else
            filename = strcat(basename,variable(k));
            [output.t,output.signal] = eprload(convertStringsToChars(filename));
        end
    end
    
    if iscell(output.t)
        tlen   = length(output.t{2});
        data.t = output.t{1};
        currrec(told:told+tlen-1)   = output.t{2}/1000;
    else
        tlen   = 1;
        data.t = output.t;
        if trec(k) ~= 0
            currrec(told:told+tlen-1) = trec(k);
        end
    end
    currsig(:,told:told+tlen-1) = output.signal;
    told = told + tlen;
end
data.signal = currsig;
data.trec   = currrec;
data.t      = data.t/1000;
end