function full_filename = create_aligned_training_dataset(varargin)
% create_aligned_training_dataset  Align a training fUS dataset to a new
% neurovascular map.
%
% INPUTS:
%   varargin
%       buffer_windows:     scalar double; How many frames in buffer?
%       n_timepoints:       scalar double; # of images fed to classifer
%       session_run_list:   (n x 2) double; [session run]
%      
% OUTPUTS:
%   full_filename:          string; path and name of the aligned dataset.

%% Variable inputs
p = inputParser;
p.addParameter('buffer_windows', 60);
p.addParameter('n_timepoints', 3);
p.addParameter('session_run_list', []);
p.parse(varargin{:});
inputs = p.Results;


%% User settings

% In the fUS-BMI, how big is the buffer window?
buffer_windows = inputs.buffer_windows; % # of images in buffer

% In the fUS-BMI, how many timepoints do we feed into the classifier?
n_timepoints = inputs.n_timepoints; % # of images fed to classifier

%% Load the data and current session neurovascular images
if isempty(inputs.session_run_list)
    gui_title = 'Select first dataset to be aligned.';
else
    gui_title = sprintf('Please choose dataset which we want to align to the session already loaded, S%dR%d', inputs.session_run_list(1), inputs.session_run_list(2));
end
data_old = load_doppler_data('gui_title', gui_title);
neurovascular_map_old = data_old.neurovascular_map;
session_run_list_old = data_old.session_run_list;

% Load the neurovascular map from today's RT recording session
data_new = load_doppler_data('session_run_list', inputs.session_run_list, ...
    'variables_to_load', {'neurovascular_map', ...
    'session_run_list'}, ...
    'gui_title', 'Select new session to align previously selected data to:');
neurovascular_map_new = data_new.neurovascular_map;
session_run_list_new = data_new.session_run_list;

% If the two images are different heights, then crop to be the same size.
% Since we always used 128 elements of the transducer, they will always be
% the same width.
if any(size(neurovascular_map_new) ~= size(neurovascular_map_old))
    % Resize the two images
    max_height = min(size(neurovascular_map_new, 1), size(neurovascular_map_old, 1));
    max_width = min(size(neurovascular_map_new, 2), size(neurovascular_map_new, 2));
    
    neurovascular_map_new_resized = neurovascular_map_new(1:max_height, 1:max_width);
    neurovascular_map_old_resized = neurovascular_map_old(1:max_height, 1:max_width);
else
    neurovascular_map_new_resized = neurovascular_map_new;
    neurovascular_map_old_resized = neurovascular_map_old;
end

% Align old neurovascular map to new neurovascular map
alignment_old_to_new = align_neurovascular_maps(neurovascular_map_new_resized, neurovascular_map_old_resized, ...
     'neurovascular_map1_session_run', session_run_list_new,...
     'neurovascular_map2_session_run', session_run_list_old);

% Convert to MATLAB format for the alignment transform
tform = affine2d(alignment_old_to_new);

% Align image to match the neurovascular map from the current recording session.
aligned_dop = imwarp(...
    data_old.dop, ...
    tform, ...
    'OutputView', imref2d(size(neurovascular_map_new_resized)));

%% Creating the training record from the previous session
% We use the previous n-1:0 timepoints before a prediction appears, where
% `n` is number of timepoints used in prediction. There are
% more advanced ways to handle this, but YAGNI for now.

prediction_ind = find(~isnan(data_old.predicted_labels) & data_old.predicted_labels ~= 0);

% Initialize empty train and train_labels variables
train = [];
train_labels = data_old.actual_labels(prediction_ind);
pb = waitbar(0, 'processing data');
for i = 1:length(prediction_ind)
    t_ind = false(size(data_old.time_acq));
    start_ind = max(1, prediction_ind(i)-buffer_windows);
    t_ind(start_ind:prediction_ind(i)) = true;
    
    % Preprocess data_old
    data_buff = aligned_dop(:, :, t_ind);
    data_buff = preprocess_data(data_buff,'zscore',true);
    [m, n, ~] = size(data_buff);
    
    % Store the last n timepoints before the prediction (including the
    % prediction timepoint)
    new_train = reshape(data_buff(:, :, end-n_timepoints+1:end), m, n*3);
    train = cat(1, train, reshape(new_train, 1, []));
    waitbar(i/length(prediction_ind));
end
close(pb);

%% Save aligned data to file

% Save as single precision to half the size. This helps initialize the
% fUS-BMI faster.
train_single = single(train);
train_labels_single = single(train_labels);

% Specify where you want to save the data
title_string = 'Where do you want to save the full aligned data?';
disp(title_string);
savepath = uigetdir(get_data_path('path_type', 'aligned'), title_string);

full_filename = fullfile(savepath, ...
    sprintf('%s_aligned_S%dR%d_to_S%dR%d_training_data_full.mat', ...
        datestr(now(), 'yymmddHHMMSS'), ...
        session_run_list_old(1), session_run_list_old(2), ...
        session_run_list_new(1), session_run_list_new(2)));

% Save the data to the desired location
save(full_filename, ...
    'train_single', 'train_labels_single', 'session_run_list_old', 'session_run_list_new', ...
    'tform', 'neurovascular_map_old', 'neurovascular_map_new');

end