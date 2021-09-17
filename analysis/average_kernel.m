function [] = average_kernel(input,output)

start = strfind(input,'CH3');
stop  = start + length('CH3a');
part1 = input(1:start-1);
part2 = input(stop:end);

CH3groups = ["CH3a","CH3b","CH3c","CH3d"];

if isempty(start)
    currK   = load(input);
    Ka   = currK.K.td;
    Ktime   = Ka;
else
    inputa = strcat(part1,CH3groups(1),part2);
    inputb = strcat(part1,CH3groups(2),part2);
    inputc = strcat(part1,CH3groups(3),part2);
    inputd = strcat(part1,CH3groups(4),part2);
    currK = load(inputa);
    currKb = load(inputb);
    currKc = load(inputc);
    currKd = load(inputd);
    Ka = currK.Ka;
    Kb = currKb.Kb;
    Kc = currKc.Kc;
    Kd = currKd.Kd;
    Ktime = (Ka + Kb + Kc + Kd)/4;
end
vt   = currK.vt;
V3   = currK.V3;
T    = currK.T;
Sys  = currK.Sys;
Exp  = currK.Exp;
Opt  = currK.Opt;

kernel = Ktime;

save(output,'kernel','vt','V3','T','Sys','Exp','Opt');


end