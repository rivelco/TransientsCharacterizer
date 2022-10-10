% PAPS orderer - This function takes as arguments the index of an esamble
% and the PAPS_INDEXED matrix, where all the PAPS scores from each ensamble
% are stored. The function returns the scores from a given ensamble with
% index ensIdx and returns the idx associated to that arrange.
function [scores, idx] = PAPSorderer(ensIdx, PAPS_INDEXED)
    % Here we retrive the data from the ensIdx ensamble
    ensData = PAPS_INDEXED(:, ensIdx);
    % Now convert the data to a matrix, in the first row we'll find the
    % indexes
    idx = cell2mat(ensData(1,1));
    % In the second row we'll find the scores
    scores = cell2mat(ensData(2,1));
    % Now sort by using the scores, from max to min, also conservating the
    % index from each element (in I variable) to sort the indexes
    [scores, I] = sort(scores, 'descend');
    % This sorts the index array while keeping the order from scores
    idx = idx(I);
end

