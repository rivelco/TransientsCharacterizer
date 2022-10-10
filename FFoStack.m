function [sP, mP, eP] = FFoStack(activity, spikes, window, doPlot)
    
    [dataNorm, startPoint, maxPoint, ~] = FFoTimePoints(activity, spikes);

    addWindow = zeros(window, 1);
    activity = [addWindow; dataNorm; addWindow];

    if doPlot
        figure(1)
    end
    
    snaps = zeros(window*2, length(maxPoint));
    for i = 1:length(maxPoint)
        mpNow = maxPoint(i)+window;
    
        currSp = mpNow - window;
        currEp = mpNow + window - 1;
        
        snap = activity(currSp:currEp);
        snaps(:, i) = snap;
        
        if doPlot
            %plot(snap, 'Color', [0.7 0.7 0.7]);
            plot(snap)
            hold on
        end
    end
    
    tauUp = maxPoint - startPoint;
    
    strFin = window-round(mean(tauUp));
    
    spike = zeros(2*window, 1);
    spike(strFin:window) = 1;
    
    final = mean(snaps, 2);
    
    [~, sP, mP, eP] = FFoTimePoints(final, spike);

    if doPlot
        plot(final, 'color', 'blue', 'LineWidth', 3)
        xline(strFin, 'Color', 'green', 'LineWidth', 2)
        xline(window, 'Color', 'blue', 'LineWidth', 2)
        xline(mean(eP), 'Color', 'red', 'LineWidth', 2)
    end
end
    
    

