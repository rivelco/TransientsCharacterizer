function [dataNorm, startPoint, maxPoint] = FFoTimePoints(FFo, spikes)
    % Get the start and end points of spikes, the endpoint is also the top
    % this is given by Aratea software
    startPoint = getStartpoints(spikes);
    maxPoint = getEndpoints(spikes);
    
    % Normalize the data, using the maximum value and divide everything by
    % that
    FFo = FFo - min(FFo);
    maxFFo = max(FFo);
    dataNorm = FFo/maxFFo;

    if isempty(startPoint)
        startPoint = NaN;
        maxPoint = NaN;
        return
    end
end