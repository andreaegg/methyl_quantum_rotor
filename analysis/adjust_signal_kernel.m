function [t,tK,signal,K] = adjust_signal_kernel(t,tK,signal,K)

factor = mean(diff(tK)) / mean(diff(t));
lK = size(K,2);
lS = length(signal);

if lK < lS
    currt = t(1:lK);
    if t(1) < tK(1)
        for n = 1:lK
            start  = find(round(currt,3) == round(tK(n),3));
            if ~isempty(start)
                break
            end
        end
        t      = t(start:end);
        signal = signal(start:end);
        tK     = tK(n:end);
        K      = K(:,n:end);
    elseif t(1) > tK(1)
        for n = 1:lK
            start = find(round(tK,3) == round(currt(n),3));
            if ~isempty(start)
                break
            end
        end
        tK     = tK(start:end);
        K      = K(:,start:end);
        t      = t(n:end);
        signal = signal(n:end);
    end
elseif lK > lS
    currK = tK(1:lS);
    if t(1) < tK(1)
        for n = 1:lS
            start  = find(round(t,3) == round(currK(n),3));
            if ~isempty(start)
                break
            end
        end
        t      = t(start:end);
        signal = signal(start:end);
        tK     = tK(n:end);
        K      = K(:,n:end);
    elseif t(1) > tK(1)
        for n = 1:lS
            start = find(round(currK,3) == round(t(n),3));
            if ~isempty(start)
                break
            end
        end
        tK     = tK(start:end);
        K      = K(:,start:end);
        t      = t(n:end);
        signal = signal(n:end);
    end
else
    if t(1) < tK(1)
        start  = find(round(t,3) == round(tK,3));
        start  = start(1);
        t      = t(start:end);
        signal = signal(start:end);
        tK     = tK(start:end);
        K      = K(:,start:end);
    elseif t(1) > tK(1)
        start  = find(round(tK,3) == round(t,3));
        start  = start(1);
        t      = t(start:end);
        signal = signal(start:end);
        tK     = tK(start:end);
        K      = K(:,start:end);
    end
end



if factor == 1
    if length(t) > length(tK)
        ende = find(round(t,3) == round(tK(end),3));
        t      = t(1:ende);
        signal = signal(1:ende);
    elseif length(t) < length(tK)
        ende = find(round(tK,3) == round(t(end),3));
        tK   = tK(1:ende);
        K    = K(:,1:ende);
    end
elseif factor > 1
    t      = t(1:factor:end);
    signal = signal(1:factor:end);
    if length(t) > length(tK)
        ende = find(round(t,3) == round(tK(end),3));
        t      = t(1:ende);
        signal = signal(1:ende);
    elseif length(t) < length(tK)
        ende = find(round(tK,3) == round(t(end),3));
        tK   = tK(1:ende);
        K    = K(:,1:ende);
    end    
elseif factor < 1
    tK = tK(1:1/factor:end);
    K  = K(:,1:1/factor:end);
    if length(t) > length(tK)
        ende = find(round(t,3) == round(tK(end),3));
        t      = t(1:ende);
        signal = signal(1:ende);
    elseif length(t) < length(tK)
        ende = find(round(tK,3) == round(t(end),3));
        tK   = tK(1:ende);
        K    = K(:,1:ende);
    end
    
else
    error('Not signal but kernel must be cropped / adjusted');
end

end