clear, clc

% CRFs file and K for significant vectors
CRFsFile = "Stoixeion_01_CRFS";
K = "3";

% AUC and Node Strength file
AUCnNSfile = "Stoixeion_01";

% FFo file
load("ffo\Stoixeion_01.mat")

% Now load the results from the CRF analysis, first from the phi 11
load("crf/" + CRFsFile + "_K" + K + " 11/results.mat", 'PAPS_INDEXED');
PAPS11 = PAPS_INDEXED;
% And later for the phi 10
load("crf/" + CRFsFile + "_K" + K + " 10/results.mat", 'PAPS_INDEXED');
PAPS10 = PAPS_INDEXED;

% Now load the results from the CRF analysis, AUCs ans node strength, 
% first from the phi 11
load("aucns/CRF_" + AUCnNSfile + "_phi11.mat");
property11 = AUCs;
% And later for the phi 10
load("aucns/CRF_" + AUCnNSfile + "_phi10.mat");
property10 = nodes_strength;

% Get the number of ensembles and the number of cells for each one
numEnsembles = max(sec_Pk_frames);
cellsPerEnsemble = cell(1, numEnsembles);
scores10 = cell(1, numEnsembles);
scores11 = cell(1, numEnsembles);
idxs10 = cell(1, numEnsembles);
idxs11 = cell(1, numEnsembles);

% For every ensamble
for ensemble = 1:numEnsembles
    % Get the ID of every cell in that ensemble
    cellsID = Pools_coords(:,3,ensemble);
    % Keep only the ones that have data (the actual cells in ensamble)
    cellsID = cellsID(cellsID > 0);

    % Now sort the neurons identificator or index (idx1X) given the score
    % obtained in the CRF analysis. This is done for both the phi 10 and 11
    [score10, idx10] = PAPSorderer(ensemble, PAPS10);
    [score11, idx11] = PAPSorderer(ensemble, PAPS11);

%     [score10, idx10] = AUCsNSorderer(ensemble, property10);
%     [score11, idx11] = AUCsNSorderer(ensemble, property11);
    
    % Save the cells in that ensemble
    cellsPerEnsemble{ensemble} = cellsID;
    % Save the score of the cells identified with CRFs with phi11 and phi10
    scores10{ensemble} = score10; 
    scores11{ensemble} = score11;
    % Save the ID of every cell identified by CRFs with phi11 and phi10
    idxs10{ensemble} = idx10; 
    idxs11{ensemble} = idx11; 
end

% Get the number of frames with activity
framesActiv = size(Pks_Frame, 2);
% Get the number of frames registered
frames = size(Spikes, 2);

% This is another way to store the activity in time of every ensamble 
actPerEnsemble = zeros(numEnsembles, frames);

%Para hacer la grafica de estados en funcion del numero de frames 
figure(3); clf; 
hold on
set(gcf,'color','w');
for hh = 1:framesActiv
    % True if this activity is from any ensamble
    if sec_Pk_frames(hh) > 0
        % Mark that this ensemble had activity in that frame
        actPerEnsemble(sec_Pk_frames(hh), Pks_Frame(hh)) = 1;
        % Plot with an empty blue circle
        plot(Pks_Frame(hh),sec_Pk_frames(hh), 'ob')
    end
end
% Set the limits, only for display purposes
xlim([1 frames]); ylim([0.5 max(sec_Pk_frames)+0.5])
set(gca,'ytick',1:framesActiv);
box on
xlabel('frame'); ylabel('ensemble')

ensembleQuery = 5;
PSNs = idxs10{ensembleQuery};
I = intersect(cellsPerEnsemble{ensembleQuery}, PSNs);

numPSN = length(PSNs);
for cell = 1:1
    currPSN = PSNs(cell);
    itsSpikes10 = Spikes(currPSN, :);
    for frame = 1:length(itsSpikes10)
        if itsSpikes10(frame)
            xline(frame, 'Color', 'red')
        end
    end
end

figure(4), clf
% Iterate over all the ensambles
for ensemble = 1:numEnsembles
    % Get the PSNs and PCNs for that ensemble and the number of each one
    PSNs = idxs10{ensemble};
    PCNs = idxs11{ensemble};
    numPSN = length(PSNs);
    numPCN = length(PCNs);
    
    % To compare spike activity vs ensemble activity
    andAC10 = zeros(1, frames);
    contAC10 = zeros(1, numPSN);
    contWA10 = zeros(1, numPSN);
    propAC10 = zeros(1, numPSN);
    orSpikes10 = zeros(1, frames);
    
    andAC11 = zeros(1, frames);
    contAC11 = zeros(1, numPCN);
    contWA11 = zeros(1, numPCN);
    propAC11 = zeros(1, numPCN);
    orSpikes11 = zeros(1, frames);
    
    % Get the binary activity vector for the current ensemble
    actEnsemble = actPerEnsemble(ensemble, :);
    % Iterate over all of its PSNs (and PCNs)
    for cell = 1:numPSN
        % Get the current cell for each set
        currPSN = PSNs(cell);
        currPCN = PCNs(cell);
        
        % Get the spikes for the PSN
        itsSpikes10 = Spikes(currPSN, :);
        % Load that activity to the orSpikes using or (only adds activity)
        orSpikes10 = bitor(itsSpikes10, orSpikes10);
        % Now create the binary vector of matches Spikes vs Ensemble
        % activity using and
        andAC10 = bitand(orSpikes10, actEnsemble);
        % The total number of matches will be the sum of that vector
        contAC10(cell) = sum(andAC10);
        % The total nomber of mistakes will be total number of frames with
        % spike activity minus the number of correct matches
        contWA10(cell) = sum(orSpikes10) - contAC10(cell);
        % Calculate the proportion of correct/wrong matches
        propAC10(cell) = contAC10(cell)/contWA10(cell);
    
        % Perform the same exact thing with the PCNs
        itsSpikes11 = Spikes(currPCN, :);
        orSpikes11 = bitor(itsSpikes11, orSpikes11);
        andAC11 = bitand(orSpikes11, actEnsemble);
        contAC11(cell) = sum(andAC11);
        contWA11(cell) = sum(orSpikes11) - contAC11(cell);
        propAC11(cell) = contAC11(cell)/contWA11(cell);

        % With every iteration the cells spikes vector will be containing
        % the activity of all the neurons in every set (PSNs or PCNs)
    end

    subplot(numEnsembles, 1, numEnsembles-ensemble+1)
    plot(1:numPSN, propAC10, 'Color', 'red')
    hold on
    plot(1:numPSN, propAC11, 'Color', 'blue')
end


% for ensemble = 1:numEnsembles
%     a = idxs10{ensemble};
%     b = idxs11{ensemble};
%     c = intersect(a, b);
%     disp("phi 10: " + length(a) + "; phi 11: " + length(b) + ...
%         "; Intersect: " + length(c) );
% end
% 
% for ensemble = 1:numEnsembles
%     a = idxs10{ensemble};
%     b = cellsPerEnsemble{ensemble};
%     c = intersect(a, b);
%     disp("phi 10: " + length(a) + "; Ens: " + length(b) + ...
%         "; Intersect: " + length(c) );
% end
% 
% for ensemble = 1:numEnsembles
%     a = idxs11{ensemble};
%     b = cellsPerEnsemble{ensemble};
%     c = intersect(a, b);
%     disp("phi 11: " + length(a) + "; Ens: " + length(b) + ...
%         "; Intersect: " + length(c) );
% end

