function [dataNorm, startPoint, maxPoint, endPoint] = FFoTimePoints(activity, spikes)
    % Retrieve the data 
    data = activity;

    % Get the start and end points of spikes, the endpoint is also the top
    % this is given by Aratea software
    startPoint = getStartpoints(spikes);
    maxPoint = getEndpoints(spikes);
    
    % Now we retrieve the size of the data
    [rows, ~] = size(data);
    
    % Normalize the data, using the maximum value and divide everything by
    % that
    maxFFo = max(data);
    dataNorm = data/maxFFo;

    % Calculate the basal mean, using the mean of activity before the first
    % spike
    %basal = mean(dataNorm(1:startPoint(1)));
    basal = dataNorm(startPoint(1));
    
    % Here we calculate the timepoints where the spike starts and decreases
    % 63% so we can calculate lambda
    % Get the maximum value given by EP
    maxSpk = zeros(length(maxPoint), 1);
    neededSpk = zeros(length(maxPoint), 1);
    endPoint = zeros(1, length(maxPoint));
    % Iterate over all the spikes identified with Aratea
    for i = 1:length(maxPoint)
        % Get the current temporal index, now we're standing in the top of
        % the AP
        currIdx = maxPoint(i);
        % Get the value of that action potencial, right at the top
        currVal = dataNorm(currIdx);
        % Register this value as the max of the AP
        maxSpk(i) = currVal;
        % Now we calculate "how much spike" is needed to be considered that
        % the activity decreased 63%. So we fubstract the basal to our AP
        % max, then multiply by 0.37 to get the value where the activity
        % decreased 63%. We also sum the basal at the end.
        neededSpk(i) = ((maxSpk(i)-basal)*.37)+basal;
        % Now we iterate over all the activity looking for the point where
        % the activity is less than the calculated before
        while currVal > neededSpk(i) && currIdx < rows
            % If we're still above that value we prepare to check the next
            % one, updating the index and also the value to check
            currIdx = currIdx + 1;
            currVal = dataNorm(currIdx);
        end
        % When we're out of that loop, means that we found that time spot
        % or that we reached the end of file
        endPoint(i) = currIdx;
    end
    
    % Plot with lines
%     plot(dataNorm);
%     xline(startPoint, '-g')
%     xline(maxPoint, '-r')
%     xline(endPoint, '-b')
%     yline(basal)

%     nBins = rows/bin;
%     dataBins = zeros(nBins, 1);
%     for i = 1:nBins
%         ini = ((i-1)*bin)+1;
%         fin = i*bin;
%         dataBins(i, 1) = mean(dataNorm(ini:fin, 1));
%     end
%     
%     dataDif = zeros(nBins, 1);
%     for i = 1:nBins
%         if i == 1
%             prev = 1;
%         else
%             prev = i-1;
%         end
%         dataDif(i) = dataBins(i) - dataBins(prev);
%     end
end