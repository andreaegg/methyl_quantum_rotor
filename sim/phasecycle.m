function [pulse_pc] = phasecycle(frequency,phase,operator)

sx = operator{1};
sy = operator{2};

if isempty(phase)
    phase = "+x";
end


for k = 1:length(phase)
    
    switch phase(k)
        case "+x"
            pulse_pc{k} = frequency*sx;
        case "-x"
            pulse_pc{k} = -frequency*sx;
        case "+y"
            pulse_pc{k} = frequency*sy;
        case "-y"
            pulse_pc{k} = -frequency*sy;
            
    end
end
end