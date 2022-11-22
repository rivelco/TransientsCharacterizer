function [tausUp, tausDown, rsqUp, rsqDown] = exponentialFit(fileName, ...
                                    threshold, Fs, spkUp, spkDown, window)
    
    disp("Reading file: " + fileName)
    load(fileName, "FFo", "Spikes");
    % window = 100;
    % Fs = 4;
    % threshold = 0.8;
    % spkUp = 50:101;
    % spkDown = 101:200;
    
    [~, cells] = size(FFo);
    disp("Total number of cells: " + cells)
    
    disp("Generating FFo stacks... ")
    cellStacks = zeros(cells, window*2);
    spikeStacks = zeros(cells, window*2);
    for cell = 1:cells
        activity = FFo(:, cell);
        spike = Spikes(cell, :);
        [cellStack, spikeStack] = FFoStack(activity, spike, window, 0);
        cellStacks(cell, :) = cellStack';
        spikeStacks(cell, :) = spikeStack';
    end
    
    modelsUp = {};
    modelsDown = {};
    tausUp = zeros(1, cells);
    tausDown = zeros(1, cells);
    rsqUp = zeros(1, cells);
    rsqDown = zeros(1, cells);
    
    goodTausUp = {};
    goodTausDown = {};
    disp('Calculating taus from stacks...')
    for cell = 1:cells
        transient = cellStacks(cell, :);
        [fitUp, gofUp, fitDown, gofDown] = expFit(transient, spkUp, spkDown);
        
        tauUp = ((1/(fitUp.b))*(1/Fs))*1000;
        rsquareUp = gofUp.rsquare; 
        tauDown = ((1/(-fitDown.b))*(1/Fs))*1000;
        rsquareDown = gofDown.rsquare; 

%         [~] = plotModelData(transient, fitUp, fitDown, spkUp, spkDown);
%         pause
        
        modelsUp{cell} = fitUp;
        modelsDown{cell} = fitDown;
        tausUp(cell) = abs(tauUp)/1000;
        tausDown(cell) = abs(tauDown)/1000;
        rsqUp(cell) = rsquareUp;
        rsqDown(cell) = rsquareDown;
    
        if rsquareUp >= threshold && rsquareDown >= threshold ...
            && tausUp(cell) < 50/Fs && tausDown(cell) < 100/Fs
            goodTausUp = [goodTausUp, tausUp(cell)];
            goodTausDown = [goodTausDown, tausDown(cell)];
        end
        
    end
    
    goodTausUp = cell2mat(goodTausUp);
    goodTausDown = cell2mat(goodTausDown);
    bothTaus = [goodTausUp; goodTausDown];
    proms = mean(bothTaus);
    
    disp("Good cells (rsqrt > " + threshold + "): " + length(goodTausDown))
    
    tauScore = (tausUp+tausDown)/2;
    % idxs = 1:length(tausUp);
    % tauPerCell = [idxs', tauScore'];
    % tauPerCell = sortrows(tauPerCell, [2, 1]);
    % save("taus.mat", "tauPerCell");
    
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
