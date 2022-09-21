clear, clc, clf
load("m21_d2.mat", "Coord_active", "FFo", "Spikes");

cell = 3;

data = FFo(:, cell);
prom2 = mean(data);
sp = getStartpoints(Spikes(cell, :));
ep = getEndpoints(Spikes(cell, :));
[rows, cols] = size(data);

basal = mean(data(1:sp(1)));

figure(1)
subplot(5, 1, 1);
plot(data);
xline(sp)
xline(ep)
yline(basal)
bin = 10;

maxFFo = max(data);

dataNorm = data/maxFFo;
basal = mean(dataNorm(1:sp(1)));

% Get the maximum value given by EP
maxSpk = zeros(length(ep), 1);
neededSpk = zeros(length(ep), 1);
pointSpk = zeros(length(ep), 1);
for i = 1:length(ep)
    currIdx = ep(i);
    currVal = dataNorm(currIdx);
    maxSpk(i) = currVal;
    neededSpk(i) = ((maxSpk(i)-basal)*.37)+basal;
    while currVal > neededSpk(i) && currIdx < 1000
        currIdx = currIdx + 1;
        currVal = dataNorm(currIdx);
    end
    pointSpk(i) = currIdx;
end

disp(ep)
disp(pointSpk)

subplot(5, 1, 2);
plot(dataNorm);
xline(sp, '-g')
xline(ep, '-r')
xline(pointSpk, '-b')
yline(basal)

% Calculate the tau for downside
tauDown = abs(ep - pointSpk');
disp(tauDown);


prom = mean(dataNorm);

dataProm = dataNorm / 10;

subplot(5, 1, 3);
plot(dataProm);

nBins = rows/bin;
dataBins = zeros(nBins, 1);
for i = 1:nBins
    ini = ((i-1)*bin)+1;
    fin = i*bin;
    dataBins(i, 1) = mean(dataNorm(ini:fin, 1));
end

subplot(5, 1, 4);
plot(dataBins);
xline(sp./bin)


subplot(5, 1, 5);
histogram(dataBins);

dataDif = zeros(nBins, 1);
for i = 1:nBins
    if i == 1
        prev = 1;
    else
        prev = i-1;
    end
    dataDif(i) = dataBins(i) - dataBins(prev);
end

figure(2)
plot(dataDif);
