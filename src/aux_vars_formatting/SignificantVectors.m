function finalFile =  SignificantVectors(fileToRead, K)
    if nargin == 0
        % Change this to select a different file
        fileToRead = uigetfile('*.*','Choose file to filter');
    end

    % Load that file
    load(fileToRead);
    
    % Get the sizes from data and UDF matrix
    [~, colsD] = size(data);
    [~, colsU] = size(UDF);
    
    % Select only those population vectors that has K or more active
    % neurons
    %K = 5;
    indexes = find(sum(data, 2) >= K);
    
    % Create the new variables
    newData = zeros(length(indexes), colsD);
    newUDF = zeros(length(indexes), colsU);
    
    % Fill the new variables using only the indexes from the valid vectors
    for i = 1:length(indexes)
        newData(i, :) = data(indexes(i), :);
        newUDF(i, :) = UDF(indexes(i), :);
    end

    % Make sure that coords has 3 dimensions, so it works right away with 
    % the CRF program
    [rowsC, colsC] = size(coords);
    if colsC < 3
        coords = [coords zeros(rowsC, 1)];
    end
    
    % Change the variables, only for the name
    data = newData;
    UDF = newUDF;
    
    % Change the filename to add the K used
    filename = filename + "_K" + K;

    % Create the new file name
    [~, finalFile, ~] = fileparts(fileToRead);
    finalFile = finalFile + "_K" + K;
    
    % Save the variables in a new file, ready to use
    save(finalFile, "filename", "coords", "data", "UDF");
end