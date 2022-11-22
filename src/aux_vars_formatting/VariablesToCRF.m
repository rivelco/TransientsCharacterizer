function [data, UDF, coords, filename] = VariablesToCRF(outputfile, ...
                                            Spikes, UDFs, Coord_active)
    disp("Ordering variables for CRF file")
    
    data = Spikes;
    [frames, cells] = size(data);
    if cells > frames
        disp("Detected more cells than frames")
        disp("Assumed inverted matrix orientation, fixing that ...")
        data = data';
        [frames, cells] = size(data);
    end
    UDF = UDFs;
    [framesUDFs, ~] = size(UDF);
    if frames ~= framesUDFs
        disp("The frames in UDFs is different from frames in spikes")
        disp("Nothing saved")
    end
    coords = Coord_active;
    [cellsCords, dims] = size(coords);
    if cells ~= cellsCords
        disp("The number of cells in spikes is diifferent from coords")
        disp("Nothing saved")
    end
    if dims == 2
        disp("There were only 2 dimensions, added a new one, rand values")
        z = randi(300, cellsCords, 1);
        coords = [coords, z];
    end
    filename = outputfile;
    
%     disp("Saving file")
%     save(filename, "data", "UDF", "coords", "filename");
end