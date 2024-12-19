function [] = pmultiply_kernel(startCH3,parts,input,output)

start = strfind(input,'CH3');
stop  = start + length('CH3a');
part1 = input(1:start-1);
part2 = input(stop:end);

CH3groups = ["CH3a","CH3b","CH3c","CH3d"];
idx       = find(CH3groups == startCH3);
groups    = CH3groups(idx:(idx+parts-1));
Ktime     = [];

if isempty(start)
    currK   = load(input);
    Ka   = currK.K.td;
    Ktime   = Ka;
else
    for k = 1:length(groups)
        input = strcat(part1,groups(k),part2);
        currK = load(input);
        if isfield(currK,'Ka')
            eval('K = currK.Ka;');
        elseif isfield(currK,'Kb')
            eval('K = currK.Kb;');
        elseif isfield(currK,'Kc')
            eval('K = currK.Kc;');
        elseif isfield(currK,'Kd')
            eval('K = currK.Kd;');
        else
            eval('K = currK.K;');
        end
        if k == 1
            Ktime = ones(length(currK.V3),length(currK.T));
        end
        Ktime = Ktime.*K;
    end
end
Ktime = Ktime/max(max(Ktime));

vt   = currK.vt;
V3   = currK.V3;
T    = currK.T;
Sys  = currK.Sys;
Exp  = currK.Exp;
Opt  = currK.Opt;

kernel = Ktime;

save(output,'kernel','vt','V3','T','Sys','Exp','Opt');
end