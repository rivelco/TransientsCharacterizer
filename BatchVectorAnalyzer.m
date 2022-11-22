% Batch vectors analyzer

% Clear the environment, console and figures
clear, clc, clf

% CRFs file and K for significant vectors
CRFsFile = "Stoixeion_04_CRFS";
K = "3";

% Declare the FFo and Spikes containing file
FFoFile = "ffo/" + "Stoixeion_04.mat";
% window = 150; Fs = 4; threshold = 0.8; spkUp = 1:151; spkDown = 151:300;
% window = 150; Fs = 4; threshold = 0.8; spkUp = 100:151; spkDown = 151:250;
window = 100; Fs = 4; threshold = 0.8; spkUp = 50:101; spkDown = 101:200;
[tausUp, tausDown, rsqUp, rsqDown] = exponentialFit(FFoFile, ...
                                    threshold, Fs, spkUp, spkDown, window);
tauScore = (tausUp+tausDown)/2;
idxs = 1:length(tausUp);
tauPerCell = [idxs', tauScore'];
tauPerCell = sortrows(tauPerCell, [2, 1]);

% Bin size determines how many neurons are acumulated for every iteration
binSize = 1;
% Ensamble is the number of ensamble to analyze inside the same database,
% it's always a number
ensamble = 1;

% Load the file to analyze, this file contains the coords, UDFs and spikes
% data
load("dbs/" + CRFsFile + "_K" + K + ".mat");
% Now load the results from the CRF analysis, first from the phi 11
load("crf/" + CRFsFile + "_K" + K + " 11/results.mat", 'PAPS_INDEXED');
PAPS11 = PAPS_INDEXED;
% Ant later for the phi 10
load("crf/" + CRFsFile + "_K" + K + " 10/results.mat", 'PAPS_INDEXED');
PAPS10 = PAPS_INDEXED;

% Now we sort the neurons identificator or index (idx1X) given the score
% obtained in the CRF analysis. This is done for both the phi 10 and 11
[scores10, idx10] = PAPSorderer(ensamble, PAPS10);
[scores11, idx11] = PAPSorderer(ensamble, PAPS11);

% Temporal code to sort the idx given the tau-score from tausFromStacks
% idx10 = tauPerCell(1:length(scores10), 1);
% idx11 = tauPerCell(end-length(scores11):end, 1);

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
    % Now we get the number of cells that was
    [contA10, contN10] = VectorsAnalyzer(data, PCNs, 1);
    sA10(i) = std(contA10) / sqrt(length(contA10));
    sN10(i) = std(contN10) / sqrt(length(contN10));
    contA10 = analyzeFunction(contA10);
    contN10 = analyzeFunction(contN10);
    yA10(i) = contA10;
    xA10(i) = zeros(1, length(contA10))+i;
    yN10(i) = contN10;
    xN10(i) = zeros(1, length(contN10))+i+stp;

    [~, I] = intersect(tauPerCell, PCNs);
    newIdx = tauPerCell(I, :);
    tScore10(i) = mean(newIdx(:, 2));
    tSEM10(i) = std(newIdx(:, 2)) / sqrt(length(newIdx(:, 2)));
end

binNumber = length(idx11) - binSize + 1;

disp("Calculating activity for PCNs (phi 11)...")
yA11 = zeros(1, binNumber); xA11 = zeros(1, binNumber);
yN11 = zeros(1, binNumber); xN11 = zeros(1, binNumber);
sA11 = zeros(1, binNumber); sN11 = zeros(1, binNumber);
tScore11 = zeros(1, binNumber); tSEM11 = zeros(1, binNumber);
for i = 1:binNumber
    PCNs = idx11(1:binSize+i-1);
    [contA11, contN11] = VectorsAnalyzer(data, PCNs, 1);
    sA11(i) = std(contA11) / sqrt(length(contA11));
    sN11(i) = std(contN11) / sqrt(length(contN11));
    contA11 = analyzeFunction(contA11);
    contN11 = analyzeFunction(contN11);
    yA11(i) = contA11;
    xA11(i) = zeros(1, length(contA11))+i;
    yN11(i) = contN11;
    xN11(i) = zeros(1, length(contN11))+i+stp;

    [~, I] = intersect(tauPerCell, PCNs);
    newIdx = tauPerCell(I, :);
    tScore11(i) = mean(newIdx(:, 2));
    tSEM11(i) = std(newIdx(:, 2)) / sqrt(length(newIdx(:, 2)));
end

figure(1)
cA = [253/255 99/255 90/255];
cB = [0 57/255 92/255];

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

figure(2)
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