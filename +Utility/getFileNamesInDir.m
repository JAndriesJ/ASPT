function [ListOfFiles, nbFiles] = getFileNamesInDir(LoadDirectory, extension)
% This function returns the names of .mat files in a given directory.
% Input: directory location as string e.g. '/home/andries/DataForProjects/SleepData/Patient1/'
% Output: cell with string entries e.g.  {'pat_1_night_1_interval_1.mat'}...                                             .   
    if 7~=exist(LoadDirectory,'dir')
        error('This directory does not exist. Make sure that you are using windows or linux filing system.')
    end
    addpath(genpath(LoadDirectory))
    ListOfFiles = dir([LoadDirectory, extension]);
    ListOfFiles = {ListOfFiles.name}';
    nbFiles = length(ListOfFiles);
end