function UDFs = TimeCourseToUDFS(sec_Pk_Frame, Pks_Frame, FFo)    
    disp("Converting time courses to UDFs ... ")
    %sec_Pk_Frame = sec_Pk_frames;
    nEnsembles  = max(sec_Pk_Frame);
    [frames, ~] = size(FFo);
    
    UDFs = zeros(frames, nEnsembles);
    sec_Pk_Frame = sec_Pk_Frame';
    
    for i = 1:length(sec_Pk_Frame)
        frame = Pks_Frame(i);
        ensemble = sec_Pk_Frame(i);
        if ensemble
            UDFs(frame, ensemble) = 1;
        end
    end    
    disp("Done converting")
end
