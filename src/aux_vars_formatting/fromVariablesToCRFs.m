clear, clc

filesToAnalize = ["Stoixeion_01_04b",...
"Stoixeion_01_04c",...
"Stoixeion_01b",...
"Stoixeion_01c",...
"Stoixeion_04_05b",...
"Stoixeion_04_05c",...
"Stoixeion_04b",...
"Stoixeion_04c",...
"Stoixeion_05b",...
"Stoixeion_05c"];

directTC = 'dbs/generadores/';
directFFo = 'dbs/generadores/';
timecourseExt = '';
finalExt = '.mat';

for fn = filesToAnalize
    disp("Reading file: " + fn)

    FFoFile = strcat(directFFo, fn, finalExt);
    disp("Loading spike-data containing file: " + FFoFile)
    load(FFoFile, 'Spikes', 'Coord_active', 'FFo')
    
    disp("Generating UDFs...")
    
    timecourseFile = strcat(directTC, fn, timecourseExt, finalExt);
    disp("Loading time-course containing file: " + timecourseFile)
    load(timecourseFile, 'sec_Pk_frames', 'Pks_Frame');
    UDFs = TimeCourseToUDFS(sec_Pk_frames, Pks_Frame, FFo);
%     disp("Saving UDFs.mat file with UDFs")
%     save("UDFs", "UDFs");

    outputfile = strcat(fn, "_CRFS");
    [data, UDF, coords, filename] = VariablesToCRF(outputfile, ...
                                    Spikes, UDFs, Coord_active);
    disp("Saving file")
    save(filename, "data", "UDF", "coords", "filename");
    disp("Done saving CRFs file")
    
    disp("Saving significant vectors")
    SignificantVectors(outputfile, 0);
    SignificantVectors(outputfile, 3);
    SignificantVectors(outputfile, 4);
    disp("Done saving significant vectors")
%     pause
end