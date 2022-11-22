function [cellPlot, spike] = FFoStack(activity, spikes, window, threshold) 
    [dataNorm, startPoint, maxPoint] = FFoTimePoints(activity, spikes);

    addWindow = zeros(window, 1);
    activity = [addWindow; dataNorm; addWindow];
    
    snaps = zeros(window*2, 1);
    it = 1;
    for i = 1:length(maxPoint)
        mpNow = maxPoint(i)+window;
    
        currSp = mpNow - window;
        currEp = mpNow + window - 1;

        if isnan(mpNow)
            snap = zeros(2*window, 1)+1;
        else 
            snap = activity(currSp:currEp);
        end
        
        normSnap = snap/max(snap);
        spkUp = 50:101;
        spkDown = 101:200;

        if threshold > 0
            [~, gofUp, ~, gofDown] = expFit(normSnap', spkUp, spkDown);
            rsquareUp = gofUp.rsquare; 
            rsquareDown = gofDown.rsquare; 
            if rsquareUp >= threshold && rsquareDown >= threshold
                snaps(:, it) = normSnap;
                it = it + 1;
            end 
        else
            snaps(:, it) = normSnap;
            it = it + 1;
        end 
    end
    lengthUp = maxPoint - startPoint;
    iniSpike = window-round(mean(lengthUp));
    
    spike = zeros(2*window, 1);
    if ~isnan(iniSpike)
        spike(iniSpike:window) = 1;
    end
    
    cellPlot = mean(snaps, 2);
end
    
    

