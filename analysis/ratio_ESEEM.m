function [t,ratio,signalH,signalD] = ratio_ESEEM(factor,tH,signalH,tD,signalD,plotflag)

if nargin < 6
    plotflag = true;
end

if plotflag
    figure('Position',[200 200 400 300])
    plot(tH,real(signalH),'k-','LineWidth',1.2);
    hold on;
    plot(tD,real(signalD),'Color',[0.6 0.6 0.6],'LineWidth',1.2);
    xlabel('T [$\mu$s]','Interpreter','latex')
    ylabel('intensity [a.u.]','Interpreter','latex')
end

if factor == 1
        if length(tH) < length(tD)
            tD = tD(1:length(signalH));
            signalD = signalD(1:length(signalH));
        elseif length(tH) > length(tD)
            tH = tH(1:length(signalD));
            signalH = signalH(1:length(signalD));
        end
        half = length(tD)/2;
        stop = length(tD);
        ratio = signalH./signalD;
        t = tH;
    elseif factor > 1
        signalD = signalD(1:factor:end);
        tD = tD(1:factor:end);
        if length(tH) < length(tD)
            tD = tD(1:length(signalH));
            signalD = signalD(1:length(signalH));
            half = length(tH)/2;
            stop = length(tH);
        elseif length(tH) > length(tD)
            tH = tH(1:length(signalD));
            signalH = signalH(1:length(signalD));
            half = length(tD)/2;
            stop = length(tD);
        end
        if tH(1) == tD(1)
            ratio = signalH./signalD;
            t = tH;
        end
    elseif factor < 1
        signalH = signalH(1:1/factor:end);
        tH = tH(1:1/factor:end);
        if length(tH) < length(tD)
            tD = tD(1:length(signalH));
            signalD = signalD(1:length(signalH));
            half = length(tH)/2;
            stop = length(tH);        
        elseif length(tH) > length(tD)
            tH = tH(1:length(signalD));
            signalH = signalH(1:length(signalD));
            half = length(tD)/2;
            stop = length(tD);
        end
        if tH(1) == tD(1)
            ratio = signalH./signalD;
            t = tH;
        end
    end
    
    difference = mean(signalH(half:stop) - signalD(half:stop));
    signalD = signalD + difference;
    
    if plotflag
        plot(tD,real(signalD),'b','LineWidth',1.2)
        xlim([0 max(t)])
        legend({'exp. NO-H','exp. NO-D','corr. NO-D'},'Interpreter','latex');
        set(gca,'FontSize',13)
    end



end