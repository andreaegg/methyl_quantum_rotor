function [] = average_kernel(input,output)

start = strfind(input,'CH3');
stop  = start + length('CH3a');
part1 = input(1:start-1);
part2 = input(stop:end);

CH3groups = ["CH3a","CH3b","CH3c","CH3d"];
Ktime = 0;
Kfreq = 0;

if isempty(start)
    currK   = load(input);
    Ka   = currK.K.td;
    Ktime   = Ka;
else
    for m = 1:length(CH3groups)
        
        currK   = load(strcat(part1,CH3groups(m),part2));
        
        Ktime = Ktime + currK.K.td/length(CH3groups);

    end
end
vt   = currK.vt;
V3   = currK.V3;
T    = currK.T;
Sys  = currK.Sys;
Exp  = currK.Exp;
Opt  = currK.Opt;

K.time = Ktime;

save(output,'K','vt','V3','T','Sys','Exp','Opt');


end