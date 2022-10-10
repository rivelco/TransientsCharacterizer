function [contA,contN] = VectorsAnalyzer(data, PCNs, K)
    %clear, clc
    
    % Now we read the data table, we need tho find the vectors where at least K
    % PCNs were active
    [rows, cols] = size(data);

    % To get random data
    %PCNs = randi(cols, 4, 1);
    
    % Sort the array so the algorithm work properly
    PCNs = sort(PCNs);
    
    % We're going to store the index of those vectors that has at least K
    % active neurons at the same time from the PCNs array
    activesVector = zeros(rows, 1);
    NonAVector = zeros(rows, 1);

    contN = {};
    contA = {};
    
    %% Get index from all valid vectors
    % First we iterate over all the rows
    for i = 1:rows
        % cont is only for counting
        contPCNs = 0;
        contNonPCNs = 0;
        PCNsIt = 1;
        % For each row we iterate over all the columns
        for j = 1:cols
            % Every time a 1 is found, we check if it's from a PCN or not
            if data(i, j) == 1
                % Here we iterate over all the PNCs until a match if found or
                % until we check all the PCNs available
                while PCNs(PCNsIt) < j && PCNsIt < length(PCNs)
                    PCNsIt = PCNsIt + 1; % Iterate until a match if found
                end
                % Now check what happend there
                if j == PCNs(PCNsIt)
                    % Here a PNC was found, so we count it
                    contPCNs = contPCNs + 1;
                else
                    % Here we didn't found any match of the current colum and PCN
                    % so the  activity found is from another neuron
                    contNonPCNs = contNonPCNs + 1;
                end
                % Given the fact that all the PCNs are already sorted, if we
                % found a match, the next match is going to be only after the
                % one that we already found, so the iterator is only restored
                % after the whole row has been checked
            end
        end
        % Save the activity frequency for each type of neuron
        activesVector(i) = contPCNs;      % PCNs
        NonAVector(i) = contNonPCNs;            % Non PCNs
    
        % Here we saved the 
        if contPCNs >= K
            contA = [contA, NonAVector(i)];
        else
            contN = [contN, NonAVector(i)];
        end
    end
    
    contN = cell2mat(contN);
    contA = cell2mat(contA);
end