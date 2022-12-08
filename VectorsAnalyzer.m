% VectorsAnalyzer calculates the number of cells coactive with K or more
% cells from the PRNs vector for a whole activity matrix.
% Input
% data -> Binary activity matrix [frames, cells]
% PRNs -> Array of cells to check for coactivity (stands for Pattern
%         Related Neurons)
% K    -> Threshold for countig as coactivation or not, if K or more PRNs
%         have activity, then the rest of activity is counted as coactivity
% Output
% prnCoactives -> Array of the amount of cells that are coactives with PRNs
% prnNoncoacts -> Array of the amount of cells that are not coactives 
% Note
% At the end this arrays togheter are going to have the same number of
% frames that have the data matrix

function [prnCoactives,prnNoncoacts] = VectorsAnalyzer(data, PRNs, K)
    % Now we read the data table, we need tho find the vectors where at 
    % least K PCNs were active

    % Get the number of frames and cells in the binary activity matrix
    [frames, cells] = size(data);
    
    % Sort the array of cells to be considered so the algorithm work 
    % properly
    PRNs = sort(PRNs);
    
    % We're going to store the index of those vectors that has at least K
    % active neurons at the same time from the PCNs array
    activesVector = zeros(frames, 1);
    NonAVector = zeros(frames, 1);

    % Arrays for the amount of coactive or non coactive cells 
    prnNoncoacts = {};
    prnCoactives = {};
    
    % Get index from all valid vectors
    % First we iterate over all the population vectors
    for ppv = 1:frames
        % We're gonna be counting the number of active cells, una counter
        % is for cells in the PRNs set and the other for the cells that are
        % not in that set
        contPRNs = 0;
        contNonPRNs = 0;
        PRNsIt = 1;
        % For each population vector we iterate over all the cells
        for cell = 1:cells
            % Every time a 1 is found, we check if it's from a PRN or not
            if data(ppv, cell) == 1
                % Here we iterate over all the PRNs until a match if found
                % or until we check all the PRNs available
                while PRNs(PRNsIt) < cell && PRNsIt < length(PRNs)
                    PRNsIt = PRNsIt + 1; % Iterate until a match if found
                end
                % Now check what happend there
                if cell == PRNs(PRNsIt)
                    % Here a PRN was found, so we count it
                    contPRNs = contPRNs + 1;
                else
                    % Here we didn't found any match of the current colum
                    % and PCN so the  activity found is from another neuron
                    contNonPRNs = contNonPRNs + 1;
                end
                % Given the fact that all the PRNs are already sorted, if 
                % we found a match, the next match is going to be only 
                % after the one that we already found, so the iterator is 
                % only restored after the whole row has been checked
            end
        end
        % Save the activity frequency for each type of neuron
        activesVector(ppv) = contPRNs;      % PCNs
        NonAVector(ppv) = contNonPRNs;            % Non PCNs
    
        % Here we saved the activity, if there was at least K PCNs involved
        % then the rest of the activity is saved as coactivated, else is
        % saved as non coactive.
        if contPRNs >= K
            prnCoactives = [prnCoactives, NonAVector(ppv)];
        else
            prnNoncoacts = [prnNoncoacts, NonAVector(ppv)];
        end
    end
    % Here just convert the cell to an array
    prnNoncoacts = cell2mat(prnNoncoacts);
    prnCoactives = cell2mat(prnCoactives);
    % Handle cases were not a single coactive or not coactive cells are
    % found, so this dont return a Nan
    if isempty(prnNoncoacts)
        prnNoncoacts = 0;
    end
    if isempty(prnCoactives)
        prnCoactives = 0;
    end
end