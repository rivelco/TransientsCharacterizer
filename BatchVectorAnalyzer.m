% Batch vectors analyzer

% Clear the environment, console and figures
clear, clc

% CRFs file and K for significant vectors
CRFsFile = "Stoixeion_01_CRFS";
K = "3";

% AUC and Node Strength file
AUCnNSfile = "Stoixeion_01";

% Declare the FFo and Spikes containing file
FFoFile = "ffo/" + "Stoixeion_01.mat";
window = 100; Fs = 4; threshold = 0.8; spkUp = 50:101; spkDown = 101:200;
[tausUp, tausDown, rsqUp, rsqDown] = exponentialFit(FFoFile, ...
                                    threshold, Fs, spkUp, spkDown, window);
% tauScore = (tausUp+tausDown)/2;
tauScore = tausUp;
idxs = 1:length(tausUp);
tauPerCell = [idxs', tauScore'];
tauPerCell = sortrows(tauPerCell, [2, 1]);

% Bin size determines how many neurons are acumulated for every iteration
binSize = 1;
% Ensamble is the number of ensamble to analyze inside the same database,
% it's always a number
ensamble = 6;

% Load the file to analyze, this file contains the coords, UDFs and spikes
% data
load("dbs/" + CRFsFile + "_K" + K + ".mat");
% Now load the results from the CRF analysis, first from the phi 11
load("crf/" + CRFsFile + "_K" + K + " 11/results.mat", 'PAPS_INDEXED');
PAPS11 = PAPS_INDEXED;
% And later for the phi 10
load("crf/" + CRFsFile + "_K" + K + " 10/results.mat", 'PAPS_INDEXED');
PAPS10 = PAPS_INDEXED;

% Now load the results from the CRF analysis, first from the phi 11
load("aucns/CRF_" + AUCnNSfile + "_phi11.mat");
property11 = nodes_strength;
% And later for the phi 10
load("aucns/CRF_" + AUCnNSfile + "_phi10.mat");
property10 = nodes_strength;

% Now we sort the neurons identificator or index (idx1X) given the score
% obtained in the CRF analysis. This is done for both the phi 10 and 11
[scores10, idx10] = PAPSorderer(ensamble, PAPS10);
[scores11, idx11] = PAPSorderer(ensamble, PAPS11);

% Now we sort the neurons identificator or index (idx1X) given the score
% obtained in the CRF analysis. This is done for both the phi 10 and 11
% [scores10, idx10] = AUCsNSorderer(ensamble, property10);
% [scores11, idx11] = AUCsNSorderer(ensamble, property11);

% Temporal code to sort the idx given the tau-score from tausFromStacks
% idx10 = intersect(tauPerCell(:,1), idx10, 'stable');
% idx11 = intersect(tauPerCell(:,1), idx11, 'stable');
% idx11 = flip(idx11);
% idx10 = flip(idx10);

% Temporal code to shuffle the idx
% idx10 = idx10(randperm(length(idx10)));
% idx11 = idx11(randperm(length(idx11)));

% Point that marks a difference of at least 2 times the std from the mean
% just to know more about the data 
maxStd10 = find(scores10 > 2*std(scores10)+mean(scores10));
maxStd11 = find(scores11 > 2*std(scores11)+mean(scores11));
% Here we display a message if there's not a single score that overcomes
% the threshold
if isempty(maxStd10)
    disp("No se cumplió el criterio de corte para phi 10")
end
if isempty(maxStd11)
    disp("No se cumplió el criterio de corte para phi 11")
end

% This parameter is used in the figures, only to shift a little bit the
% points of a group plotted to the right (if value > 0) 
stp = 0;

% Here is calculated the number of bins for the phi 10 elements, given the
% number of cells found with CRF and the bin size
binNumber = length(idx10) - binSize + 1;

% Here is initialized the different arrays that are going to be used 
disp("Calculating activity for PSNs (phi 10)...")
yA10 = zeros(1, binNumber); xA10 = zeros(1, binNumber);
yN10 = zeros(1, binNumber); xN10 = zeros(1, binNumber);
sA10 = zeros(1, binNumber); sN10 = zeros(1, binNumber);
tScore10 = zeros(1, binNumber); tSEM10 = zeros(1, binNumber);
for i = 1:binNumber
    % Here is determined the PCNs set to be used now
    PCNs = idx10(1:binSize+i-1);
    % Now we get the number of cells that are coactive or not coactive with
    % PRNs
    [Coacts10, NoCoacts10] = VectorsAnalyzer(data, PCNs, 1);
    % Calculate the SEM and the result of analyze function
    sA10(i) = std(Coacts10) / sqrt(length(Coacts10));
    sN10(i) = std(NoCoacts10) / sqrt(length(NoCoacts10));
    Coacts10 = analyzeFunction(Coacts10);
    NoCoacts10 = analyzeFunction(NoCoacts10);
    
    % Now generate the points fot the scatter plot, for coactives and non
    % coactives
    yA10(i) = Coacts10;
    xA10(i) = i;
    yN10(i) = NoCoacts10;
    xN10(i) = i+stp;
    
    % This is for the tau score figure, select the tau score for the cells
    % of the PCN array, and get the mean and SEM of it 
    [~, I] = intersect(tauPerCell, PCNs);
    newIdx = tauPerCell(I, :);
    tScore10(i) = mean(newIdx(:, 2));
    tSEM10(i) = std(newIdx(:, 2)) / sqrt(length(newIdx(:, 2)));
end

% Now calculate the number of bins for the idx11 elements
binNumber = length(idx11) - binSize + 1;

disp("Calculating activity for PCNs (phi 11)...")
% Arrays to be used
yA11 = zeros(1, binNumber); xA11 = zeros(1, binNumber);
yN11 = zeros(1, binNumber); xN11 = zeros(1, binNumber);
sA11 = zeros(1, binNumber); sN11 = zeros(1, binNumber);
tScore11 = zeros(1, binNumber); tSEM11 = zeros(1, binNumber);
for i = 1:binNumber
    % Get the new PRN to be used
    PCNs = idx11(1:binSize+i-1);
    % Calculate the coactivity and no coactivity
    [Coacts11, NoCoacts11] = VectorsAnalyzer(data, PCNs, 1);
    % Calculate the SEM from the results
    sA11(i) = std(Coacts11) / sqrt(length(Coacts11));
    sN11(i) = std(NoCoacts11) / sqrt(length(NoCoacts11));
    % Result from the analyze function
    Coacts11 = analyzeFunction(Coacts11);
    NoCoacts11 = analyzeFunction(NoCoacts11);
    % Now the points for the scatter plot
    yA11(i) = Coacts11;
    xA11(i) = i;
    yN11(i) = NoCoacts11;
    xN11(i) = i+stp;
    
    % Calculate the mean and SEM for the tau score for the cells considered
    [~, I] = intersect(tauPerCell, PCNs);
    newIdx = tauPerCell(I, :);
    tScore11(i) = mean(newIdx(:, 2));
    tSEM11(i) = std(newIdx(:, 2)) / sqrt(length(newIdx(:, 2)));
end

% Create the figure
figure(1), clf
% The colors to be used, red and blue
cA = [253/255 99/255 90/255];
cB = [0 57/255 92/255];

% Results for the phi11 elements
subplot(221)
e = errorbar(xA11, yA11, sA11, 'o');
e.MarkerSize = 3;
e.MarkerFaceColor = "auto";
e.Color = cA;
e.CapSize = 3;
hold on
e = errorbar(xN11, yN11, sN11, 'o');
e.MarkerSize = 3;
e.MarkerFaceColor = "auto";
e.Color = cB;
e.CapSize = 3;
legend('Coactivas con PCNs', 'No coactivas con PCNs')
title('Actividad de neuronas con PCNs phi 11')
xlabel('Número de PCNs incluidas')
ylabel('Cantidad promedio de neuronas')

% Results for the phi10 elements
subplot(222)
e = errorbar(xA10, yA10, sA10, 'o');
e.MarkerSize = 3;
e.MarkerFaceColor = "auto";
e.Color = cA;
e.CapSize = 3;
hold on
e = errorbar(xN10, yN10, sN10, 'o');
e.MarkerSize = 3;
e.MarkerFaceColor = "auto";
e.Color = cB;
e.CapSize = 3;
legend('Coactivas con PSNs', 'No coactivas con PSNs')
title('Actividad de neuronas con PSNs phi 10')
xlabel('Número de PSNs incluidas')
ylabel('Cantidad promedio de neuronas')

% Comparision for phi11 and phi11 coactivation
subplot(223)
e = errorbar(xA11, yA11, sA11, 'o');
e.MarkerSize = 3;
e.MarkerFaceColor = "auto";
e.Color = cB;
e.CapSize = 3;
hold on
e = errorbar(xA10, yA10, sA10, 'o');
e.MarkerSize = 3;
e.MarkerFaceColor = "auto";
e.Color = cA;
e.CapSize = 3;
legend('PCNs', 'PSNs')
title('Actividad de neuronas coactivas')
xlabel('Número de neuronas incluidas')
ylabel('Cantidad promedio de neuronas')

% Results for the phi11 and phi10 no coactivation
subplot(224)
e = errorbar(xN11, yN11, sN11, 'o');
e.MarkerSize = 3;
e.MarkerFaceColor = "auto";
e.Color = cB;
e.CapSize = 3;
hold on
e = errorbar(xN10, yN10, sN10, 'o');
e.MarkerSize = 3;
e.MarkerFaceColor = "auto";
e.Color = cA;
e.CapSize = 3;
legend('PCNs', 'PSNs')
title('Actividad de neuronas no coactivas')
xlabel('Número de neuronas incluidas')
ylabel('Cantidad promedio de neuronas')

% Figure for the tau scores
figure(2), clf
e = errorbar(xA10, tScore10, tSEM10, 'o');
e.MarkerSize = 3;
e.MarkerFaceColor = "auto";
e.Color = cA;
e.CapSize = 3;
hold on
e = errorbar(xA11, tScore11, tSEM11, 'o');
e.MarkerSize = 3;
e.MarkerFaceColor = "auto";
e.Color = cB;
e.CapSize = 3;
legend('Phi 10', 'Phi 11')
title('Promedio de Tau score por categoría de neuronas incluidas')
xlabel('Número de neuronas incluidas')
ylabel('Tau score promedio (s)')


return
%% This part of the code was made for ssome scatter plots, not in use

% cA = [253/255 99/255 90/255];
% cB = [0 57/255 92/255];
% 
% figure(1)
% scatter(xA11, yA11, 10, cA, 'filled', "jitter", 'on', 'jitterAmount', 0.1);
% hold on
% scatter(xN11, yN11, 10, cB, 'filled', "jitter", 'on', 'jitterAmount', 0.1);
% legend('Coactivas con PCNs', 'No coactivas con PCNs')
% title('Actividad de células con PCNs phi 11')
% xlabel('Número de PCNs incluidas')
% ylabel('Cantidad promedio de células')
% 
% figure(2)
% scatter(xA10, yA10, 10, cA, 'filled', "jitter", 'on', 'jitterAmount', 0.1);
% hold on
% scatter(xN10, yN10, 10, cB, 'filled', "jitter", 'on', 'jitterAmount', 0.1);
% legend('Coactivas con PCNs', 'No coactivas con PCNs')
% title('Actividad de células con PCNs phi 10')
% xlabel('Número de PCNs incluidas')
% ylabel('Cantidad promedio de células')
% 
% figure(3)
% scatter(xA11, yA11, 10, cA, 'filled', "jitter", 'on', 'jitterAmount', 0.1);
% hold on
% scatter(xA10, yA10, 10, cB, 'filled', "jitter", 'on', 'jitterAmount', 0.1);
% legend('PCNs phi11', 'PCNs phi10')
% title('Actividad de células coactivas con PCNs')
% xlabel('Número de PCNs incluidas')
% ylabel('Cantidad promedio de células')
% 
% figure(4)
% scatter(xN11, yN11, 10, cA, 'filled', "jitter", 'on', 'jitterAmount', 0.1);
% hold on
% scatter(xN10, yN10, 10, cB, 'filled', "jitter", 'on', 'jitterAmount', 0.1);
% legend('PCNs phi11', 'PCNs phi10')
% title('Actividad de células no coactivas con PCNs')
% xlabel('Número de PCNs incluidas')
% ylabel('Cantidad promedio de células')