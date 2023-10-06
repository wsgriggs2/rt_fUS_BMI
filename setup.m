%% Set up code repository and data directory
% This only needs to be run once


%% 1. Add folder and subfolder to matlab search path
% Add the parent folder above the current folder and all
% subfolders of that parent folder to the MATLAB search path. Then save the
% search path.

folder = fileparts(mfilename('fullpath')); 
addpath(genpath(folder));
rmpath(genpath(fullfile(folder, '.git')));
savepath;

%% 2. Set the path to the downloaded data
path = specify_data_path;
fprintf('\nSpecified "%s" as the path to the downloaded data.\n', path);
