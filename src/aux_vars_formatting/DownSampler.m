% Downsampler script
clear, clc
fileToRead = uigetfile('*.*','Choose file to down sample');
disp("Opening file to down sample")
load(fileToRead, "FFo");

[frames, cells] = size(FFo);

window = 8;
recordSize = 2000;
disp("Converting ...")
sizeFrames = round(frames/window);
newFFo = zeros(sizeFrames, cells);
for frame = 1:window:frames
    indx = (frame-1)/window + 1;
    fin = frame + window-1;
    if fin > frames
        newFFo(indx, :) = mean(FFo(frame:end, :));
    else
        newFFo(indx, :) = mean(FFo(frame:fin, :));
    end
end
frames = size(newFFo, 1);
sp = randi(frames-recordSize);

FFoDS = newFFo(sp:sp+recordSize-1, :);

disp("Saving the new FFo file ...")
% Create the new file name
[~, finalFile, ~] = fileparts(fileToRead);
finalFile = finalFile + "_DS";
% Save the variables in a new file, ready to use
save(finalFile, "FFoDS");

disp("Generating and saving the new coords ...")
coords = randi(300, cells, 2);
save(finalFile+"_coords", "coords")

disp("Done down sampling")