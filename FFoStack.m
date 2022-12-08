% FFoStack generates the stack of the spikes from the FFo activity of a
% single cell.
% Input
% activity   -> The FFo data of the cell (vector with length of frames)
% spikes     -> Binary vector with the presence of the spike, given by
%               Aratea
% window     -> Number of frames that are going to contain the rising and
%               falling part of the stack, so the whole stack is going to
%               have twice as much frames as this window integer, the pike of
%               the spike will be in the frame [window].
% threshold  -> The minimun value of r-squared fit that a transient must
%               have to be considered as valid and to be included in the
%               stack
% Output
% cellPlot   -> Mean of every valid spike in the register, vector of
%               2*window with normalized FFo data
% spike      -> Binary matrix showing where the spike begins

function [cellPlot, spike] = FFoStack(activity, spikes, window, threshold)
    % Get the frame number of every spike start point and max point, also
    % get the normalized FFo data (FFo-min)/max
    [dataNorm, startPoint, maxPoint] = FFoTimePoints(activity, spikes);
    
    % Here we add another window at the begining and end of the activity to
    % prevent getting out of the array, that window is filled with zeros
    addWindow = zeros(window, 1);
    activity = [addWindow; dataNorm; addWindow];
    
    % This contains a window arround every spike in the register
    snaps = zeros(window*2, 1);
    it = 1;
    % Iterate over all the spikes declared in the spikes variable
    for i = 1:length(maxPoint)
        % Get the frame that keeps current max point (pike of the current 
        % spike). Here is added [window] to compensate the window added at 
        % the begining
        mpNow = maxPoint(i)+window;
        
        % Get the frame of the start point, [window] frames behind the pike
        currSp = mpNow - window;
        % Get the frame of the end point, [window] frames after the pike
        currEp = mpNow + window - 1;
    
        % This is true if there is not a single spike in the register
        if isnan(mpNow)
            % In that case the snap is just a bunch of 1's
            snap = zeros(2*window, 1)+1;
        else 
            % If the spike exists then the snap is the frames arround that
            snap = activity(currSp:currEp);
        end
        
        % Now, normalize the spike, again, so all the spikes have the same
        % heigth, so to speak
        normSnap = snap/max(snap);
        % Now set the range of frames to be fitted in the model to check if
        % that spike looks like a spike. Use the second quarter for the
        % rising part and the second half for the falling part
        spkUp = window/2:window+1;
        spkDown = window+1:window*2;
    
        % If we actually set a threshold to check, then perform the fit
        if threshold > 0
            % Get the model for the snap
            [~, gofUp, ~, gofDown] = expFit(normSnap', spkUp, spkDown);
            % Retrive the r-squared
            rsquareUp = gofUp.rsquare; 
            rsquareDown = gofDown.rsquare; 
            % If that is good enough then it is stored in snaps
            if rsquareUp >= threshold && rsquareDown >= threshold
                snaps(:, it) = normSnap;
                it = it + 1;
            end 
        else
            % Here we dont want to filter anything so the snap pass
            snaps(:, it) = normSnap;
            it = it + 1;
        end 
    end
    % Now generate the spike from that stack, use this with caution
    % Get the difference between all the max and start points
    lengthUp = maxPoint - startPoint;
    % Set the start of the spike at the mean behind the actual max [window]
    iniSpike = window-round(mean(lengthUp));
    % Generate the binary spike vector for the stack
    spike = zeros(2*window, 1);
    if ~isnan(iniSpike)
        % If there was a spike then the spike goes from the selected
        % startpoint to the window frame
        spike(iniSpike:window) = 1;
    end
    % Now calculate the mean of all the valid snaps, that will be another
    % vector, this one has the stack.
    cellPlot = mean(snaps, 2);
end
    
    

