function finalFile =  SignificantVectors(fileToRead)
    if nargin == 0
        % Change this to select a different
        fileToRead = uigetfile('*.*','Choose file to filter');
    end

    % Load that file
    load(fileToRead);
    
    % Get the sizes from data and UDF matrix
    [~, colsD] = size(data);
    [~, colsU] = size(UDF);
    
    % Select only those population vectors that has three or more active
    % neurons
    indexes = find(sum(data, 2) >= 3);
    
    % Create the new variables
    newData = zeros(length(indexes), colsD);
    newUDF = zeros(length(indexes), colsU);
    
    % Fill the new variables using only the indexes from the valid vectors
    for i = 1:length(indexes)
        newData(i, :) = data(indexes(i), :);
        newUDF(i, :) = UDF(indexes(i), :);
    end
    
    % Change the variables, only for the name
    data = newData;
    UDF = newUDF;
    
    % Create the new file name
    [~, finalFile, ~] = fileparts(fileToRead);
    finalFile = finalFile + "_sig";
    
    % Save the variables in a new file, ready to use
    save(finalFile, "filename", "coords", "data", "UDF");
end