% exponentialFit calculates the spike stack from all the cells in a
% fluorescence register and creates two exponential fits for each stack of
% each cell, one for the rising part and another for the falling part, then
% calculates the tau and r-squared for each curve.
% Input
% fileName -> The file were is contained the FFo and Spikes matrix
%          -> FFo is a matrix of [frames, cells] for the fluorescence
%          -> Spikes is a binary matrix [cells, frames] for the spikes
% threshold -> Value between 0 and 1, this marks the threshold to consider
%              an adjust as good (for r-squared from the exponential fit)
% Fs       -> Sampling frequency in Hz
% spkUp    -> Range of frames from the FFo stack to be fitted in the
%             first exponential fit
% spkDown  -> Range of frames from the FFo stack to be fitted in the
%             second exponential fit
% window   -> Number of frames that are going to contain the rising and
%             falling part of the stack, so the whole stack is going to
%             have twice as much frames as this window integer, the pike of
%             the spike will be in the frame [window].
% Output
% tausUp   -> Vector with the tau up for every cell (in seconds)
% tausDown -> Vector with the tau down for every cell (in seconds)
% rsqUp    -> r-squared value for the up fit of the stack of every cell
% rsqDown  -> r-squared value for the down fit of the stack of every cell

function [tausUp, tausDown, rsqUp, rsqDown] = exponentialFit(fileName, ...
                                    threshold, Fs, spkUp, spkDown, window)
    % Read the file
    disp("Reading file: " + fileName)
    load(fileName, "FFo", "Spikes");
    % Get the total number of cells
    [~, cells] = size(FFo);
    disp("Total number of cells: " + cells)
    % Generates the FFo stack for every cell
    disp("Generating FFo stacks... ")
    % This will contain the stack for every cell registered
    cellStacks = zeros(cells, window*2);
    % Here will be the binary spike vector of the stack
    spikeStacks = zeros(cells, window*2);
    % We'll do this for every cell
    for cell = 1:cells
        % Get the activity registered for the current cell and its spike
        activity = FFo(:, cell);
        spike = Spikes(cell, :);
        % Calculate the stack
        [cellStack, spikeStack] = FFoStack(activity, spike, window, 0);
        % Store the stack and spike of that cell
        cellStacks(cell, :) = cellStack';
        spikeStacks(cell, :) = spikeStack';
    end
    
    % Arrays to store the models, taus and r-squared, for all the cells
    modelsUp = {};
    modelsDown = {};
    tausUp = zeros(1, cells);
    tausDown = zeros(1, cells);
    rsqUp = zeros(1, cells);
    rsqDown = zeros(1, cells);
    
    % Here is stored the good cells (depends on the threshold variable)
    goodTausUp = {};
    goodTausDown = {};
    % Now calculate the exponential fit, taus and r-squared for every cell
    disp('Calculating taus from stacks...')
    for cell = 1:cells
        % Retrive the stack of the current cell
        transient = cellStacks(cell, :);
        % Calculate the exponential fit
        [fitUp, gofUp, fitDown, gofDown] = expFit(transient, spkUp, spkDown);
        
        % Calculate the tau up and down, using the b variable from the
        % model and the sampling frequency to get it on seconds. Also get
        % the r-squared
        tauUp = ((1/(fitUp.b))*(1/Fs));
        rsquareUp = gofUp.rsquare; 
        tauDown = ((1/(-fitDown.b))*(1/Fs));
        rsquareDown = gofDown.rsquare; 
        
        % Uncomment this to see the fit of every cell
        % [~] = plotModelData(transient, fitUp, fitDown, spkUp, spkDown);
        % pause
        
        % Store the model, the taus and r-squared
        modelsUp{cell} = fitUp;
        modelsDown{cell} = fitDown;
        tausUp(cell) = abs(tauUp);
        tausDown(cell) = abs(tauDown);
        rsqUp(cell) = rsquareUp;
        rsqDown(cell) = rsquareDown;
        
        % This is to save the good cells
        if rsquareUp >= threshold && rsquareDown >= threshold ...
            && tausUp(cell) < 50/Fs && tausDown(cell) < 100/Fs
            goodTausUp = [goodTausUp, tausUp(cell)];
            goodTausDown = [goodTausDown, tausDown(cell)];
        end
        
    end
    
    % Convert the good cells cell to matrix and display the number of good
    % cells
    goodTausUp = cell2mat(goodTausUp);
    goodTausDown = cell2mat(goodTausDown);
    disp("Good cells (rsqrt > " + threshold + "): " + length(goodTausDown))
    
    % This is to generate the histograms about taus
%     tauScore = (tausUp+tausDown)/2;
%     figure(1)
%     subplot(221)
%     histogram(tausUp, [0.5:0.1:2])
%     title("Tau up")
%     xlabel("tau (s)")
%     subplot(222)
%     histogram(tausDown, [0:1:30])
%     title("Tau down")
%     xlabel("tau (s)")
%     subplot(223)
%     scatter(tausUp, tausDown)
%     title("Tau up / tau down")
%     xlabel("tau up (s)")
%     ylabel("tau down (s)")
%     xlim([0.5, 2])
%     ylim([0, 30])
%     subplot(224)
%     histogram(tauScore)
%     title("tau score")
%     xlabel("mean (s)")
end
