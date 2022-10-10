clear, clc, clf

load("dbs/m21_d2.mat", "FFo", "Spikes");

[rows, cols] = size(FFo);

%tUp = zeros(cols, 1);
%tDown = zeros(cols, 1);}
tUp = {};
tDown = {};

window = 100;

cell = 0;
if cell
    activity = FFo(:, cell);
    spikes = Spikes(cell, :);
    [sP, mP, eP] = FFoStack(activity, spikes, window, true);
    return
end

for i = 1:cols
    if i == 13 || i == 43 || i == 46
        continue
    end
    cell = i;
    activity = FFo(:, cell);
    spikes = Spikes(cell, :);
    [sP, mP, eP] = FFoStack(activity, spikes, window, false);
    
    tauUp = mP - sP;
    tauDown = eP - mP;

%     tauUp = mean(tauUp);
%     tauDown = mean(tauDown);

    tUp = [tUp, tauUp];
    tDown = [tDown, tauDown];
end

tUp = cell2mat(tUp);
tDown = cell2mat(tDown);

figure(1)
c = linspace(1,10,length(tUp));
scatter(tUp, tDown, 20, c, 'filled', "jitter", 'on', 'jitterAmount', 0.3);
title('Tau promedio de cada célula')
xlabel('Tau de subida')
ylabel('Tau de bajada')

figure(2)
axis([min(tUp)-1, max(tUp)+1, min(tDown)-1, max(tDown)+1])
for k=1:length(tUp)
    text(tUp(k),tDown(k),num2str(k),'HorizontalAlignment', 'Center',...
        'VerticalAlignment', 'Middle')
end
title('Tau promedio de cada célula')
xlabel('Tau de subida')
ylabel('Tau de bajada')
