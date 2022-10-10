% Batch vectors analyzer
clear, clf
load("dbs/m21_d2_04_all_CRF_K5.mat");
load("crf/m21_d2_04_all_CRF_K5 11/results.mat", 'PAPS_INDEXED');
PAPS11 = PAPS_INDEXED;
load("crf/m21_d2_04_all_CRF_K5 10/results.mat", 'PAPS_INDEXED');
PAPS10 = PAPS_INDEXED;

% load("dbs/LCR_expert_K3.mat");
% load("crf/expertK3 11/results.mat", 'PAPS_INDEXED');
% PAPS11 = PAPS_INDEXED;
% load("crf/expertK3 10/results.mat", 'PAPS_INDEXED');
% PAPS10 = PAPS_INDEXED;

binSize = 1;
ensamble = 7;

[scores10, idx10] = PAPSorderer(ensamble, PAPS10);
[scores11, idx11] = PAPSorderer(ensamble, PAPS11);

% disp(length(scores11))
% disp(length(scores10))
% [val, pos] = intersect(scores10, scores11);
% disp(idx11')
% disp(idx10')

% Point that marks a difference of at least 2 times the std from the mean
maxStd10 = find(scores10 > 2*std(scores10)+mean(scores10));
maxStd11 = find(scores11 > 2*std(scores11)+mean(scores11));
%maxStd10 = find(scores10 > 0.75);
%maxStd11 = find(scores11 > 0.75);

if isempty(maxStd10)
    disp("No se cumplió el criterio de corte para phi 10")
end
if isempty(maxStd11)
    disp("No se cumplió el criterio de corte para phi 11")
end

%binSize = maxStd;
stp = 0;

yA10 = {}; xA10 = {};
yN10 = {}; xN10 = {};
sA10 = {}; sN10 = {};

binNumber = length(idx10) - binSize + 1;
for i = 1:binNumber    
    PCNs = idx10(1:binSize+i-1);
    [contA10, contN10] = VectorsAnalyzer(data, PCNs, 1);
    sA10 = [sA10, std(contA10) / sqrt(length(contA10))];
    sN10 = [sN10, std(contN10) / sqrt(length(contN10))];
    contA10 = analyzeFunction(contA10);
    contN10 = analyzeFunction(contN10);
    yA10 = [yA10, contA10];
    xA10 = [xA10, zeros(1, length(contA10))+i];
    yN10 = [yN10, contN10];
    xN10 = [xN10, zeros(1, length(contN10))+i+stp];
end
yA10 = cell2mat(yA10); xA10 = cell2mat(xA10);
yN10 = cell2mat(yN10); xN10 = cell2mat(xN10);
sA10 = cell2mat(sA10); sN10 = cell2mat(sN10);

yA11 = {}; xA11 = {};
yN11 = {}; xN11 = {};
sA11 = {}; sN11 = {};
binNumber = length(idx11) - binSize + 1;
for i = 1:binNumber
    PCNs = idx11(1:binSize+i-1);
    [contA11, contN11] = VectorsAnalyzer(data, PCNs, 1);
    sA11 = [sA11, std(contA11) / sqrt(length(contA11))];
    sN11 = [sN11, std(contN11) / sqrt(length(contN11))];
    contA11 = analyzeFunction(contA11);
    contN11 = analyzeFunction(contN11);
    yA11 = [yA11, contA11];
    xA11 = [xA11, zeros(1, length(contA11))+i];
    yN11 = [yN11, contN11];
    xN11 = [xN11, zeros(1, length(contN11))+i+stp];
end
yA11 = cell2mat(yA11); xA11 = cell2mat(xA11);
yN11 = cell2mat(yN11); xN11 = cell2mat(xN11);
sA11 = cell2mat(sA11); sN11 = cell2mat(sN11);

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
title('Actividad de células con PCNs phi 11')
xlabel('Número de PCNs incluidas')
ylabel('Cantidad promedio de células')

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
legend('Coactivas con PCNs', 'No coactivas con PCNs')
title('Actividad de células con PCNs phi 10')
xlabel('Número de PCNs incluidas')
ylabel('Cantidad promedio de células')

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
legend('PCNs phi11', 'PCNs phi10')
title('Actividad de células coactivas con PCNs')
xlabel('Número de PCNs incluidas')
ylabel('Cantidad promedio de células')

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
legend('PCNs phi11', 'PCNs phi10')
title('Actividad de células no coactivas con PCNs')
xlabel('Número de PCNs incluidas')
ylabel('Cantidad promedio de células')

return

figure(1)
subplot(221)
cA = [253/255 99/255 90/255];
cB = [0 57/255 92/255];
scatter(xA10, yA10, 10, cA, 'filled', "jitter", 'on', 'jitterAmount', 0.1);
hold on
scatter(xN10, yN10, 10, cB, 'filled', "jitter", 'on', 'jitterAmount', 0.1);
if ~isempty(maxStd10)
    xline(maxStd10+stp/2, 'LineWidth', 1, 'Color', [1 166/255 0]);
    legend('Coactivas con PCNs', 'No coactivas con PCNs', 'Umbral')
else
    legend('Coactivas con PCNs', 'No coactivas con PCNs')
end
title('Actividad de células con PCNs phi 10')
xlabel('Número de células incluidas')
ylabel('Células coactivas')

subplot(222)
cA = [253/255 99/255 90/255];
cB = [0 57/255 92/255];
scatter(xA11, yA11, 10, cA, 'filled', "jitter", 'on', 'jitterAmount', 0.1);
hold on
scatter(xN11, yN11, 10, cB, 'filled', "jitter", 'on', 'jitterAmount', 0.1);
if ~isempty(maxStd10)
    xline(maxStd10+stp/2, 'LineWidth', 1, 'Color', [1 166/255 0]);
    legend('Coactivas con PCNs', 'No coactivas con PCNs', 'Umbral')
else
    legend('Coactivas con PCNs', 'No coactivas con PCNs')
end
title('Actividad de células con PCNs phi 11')
xlabel('Número de células incluidas')
ylabel('Células coactivas')

subplot(223)
cA = [253/255 99/255 90/255];
cB = [0 57/255 92/255];
scatter(xA11, yA11, 10, cA, 'filled', "jitter", 'on', 'jitterAmount', 0.1);
hold on
scatter(xA10, yA10, 10, cB, 'filled', "jitter", 'on', 'jitterAmount', 0.1);
if ~isempty(maxStd10)
    xline(maxStd10+stp/2, 'LineWidth', 1, 'Color', [1 166/255 0]);
    legend('PCNs phi11', 'PCNs phi10', 'Umbral')
else
    legend('PCNs phi11', 'PCNs phi10')
end
title('Actividad de células coactivas')
xlabel('Número de células incluidas')
ylabel('Células coactivas')

subplot(224)
cA = [253/255 99/255 90/255];
cB = [0 57/255 92/255];
scatter(xN11, yN11, 10, cA, 'filled', "jitter", 'on', 'jitterAmount', 0.1);
hold on
scatter(xN10, yN10, 10, cB, 'filled', "jitter", 'on', 'jitterAmount', 0.1);
if ~isempty(maxStd10)
    xline(maxStd10+stp/2, 'LineWidth', 1, 'Color', [1 166/255 0]);
    legend('PCNs phi11', 'PCNs phi10', 'Umbral')
else
    legend('PCNs phi11', 'PCNs phi10')
end
title('Actividad de células no coactivas')
xlabel('Número de células incluidas')
ylabel('Células coactivas')