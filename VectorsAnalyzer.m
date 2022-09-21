function [contA,contN] = VectorsAnalyzer()
    %clear, clc
    
    % This array contains the index of PCNs
    % PCNs 10
    %PCNs = [2, 34, 38, 66, 67];
    %PCNs = 70;
    
    % PCNs 11
    %PCNs = [2, 3, 37, 55, 66];
    %PCNs = 55;

    load("m21_d2_04_CRF.mat");
    % PCNs for m21_d2_04_CRF 10
    %PCNs = [50, 6];
    %PCNs = [62, 67];
    %PCNs = [67, 3];
    %PCNs = [30, 65, 73];
    %PCNs = [50, 22];
    %PCNs = [67, 3, 57];
    %PCNs = [50, 6, 62, 67, 3, 30, 65, 73, 50, 22, 57];

    % PCNs for m21_d2_04_CRF 01
    %PCNs = [68, 75];
    %PCNs = [62, 69];
    %PCNs = [59, 37, 35];
    %PCNs = [65, 48];
    %PCNs = [75, 10];
    %PCNs = [58, 40];
    %PCNs = [68, 75, 62, 69, 59, 37, 35, 65, 48, 10, 58, 40];

    % PCNs for m21_d2_04_CRF 10 double
    %PCNs = [62, 45, 49];
    %PCNs = [55, 63, 73];
    %PCNs = [70, 72, 71];
    %PCNs = [57, 63, 67];
    %PCNs = [71, 74, 65];
    %PCNs = [73, 71];

    % Now we read the data table, we need tho find the vectors where at least K
    % PCNs were active
    [rows, cols] = size(data);
    
    %PCNs = randi(cols, 1, 5);
    
    PCNs = sort(PCNs);
    
    K = 2;
    
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
        contActives = 0;
        contNonA = 0;
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
                    contActives = contActives + 1;
                else
                    % He we didn't found any match of the current colum and PCN
                    % so the  activity found is from another neuron
                    contNonA = contNonA + 1;
                end
                % Given the fact that all the PCNs are already sorted, if we
                % found a match, the next match is going to be only after the
                % one that we already found, so the iterator is only restored
                % after the whole row has been checked
            end
        end
        % Save the activity frequency for each type of neuron
        activesVector(i) = contActives;      % PCNs
        NonAVector(i) = contNonA;            % Non PCNs

        if activesVector(i) >= K
            contA = [contA, NonAVector(i)];
        else
            contN = [contN, NonAVector(i)];
        end
    end
    
    contN = cell2mat(contN);
    contA = cell2mat(contA);
    
    disp(mean(contA));
    disp(mean(contN));
    
    X = [contA, contN];
    grp = [ones(size(contA)), 2.*ones(size(contN))];
    
    boxplot(X, grp)

    contN = mean(contN);
    contA = mean(contA);
end