function [train_data, train_labels] = extract_training_data_and_labels(data, varargin)
% EXTRACT_TRAINING_DATA_AND_LABELS  For each trial, find last P frames 
% preceding movement. Apply appropriate preprocessing, such as z-scoring
% and/or spatial filtering
%
% INPUTS:
%   data:                   struct; Must contain at least `behavior`,
%                           `dop`, and `timestamps`. This variable is
%                           normally loaded via LOAD_DOPPLER_DATA. Will
%                           ignore extra fields in struct.
%   varargin: 
%       zscore:             bool; Apply z-scoring?
%       spatial_filter:     (3 x 1) cell array; Specify parameters that
%                           will eventually be provided to FSPECIAL.
%                           Example: {'disk', 2, 0}.
%       training_set_size:  Scalar double; How many frames do we want to
%                           keep per trial?
%       buffer_size:        Scalar double; How many frames should be in the
%                           rolling buffer?
%
% OUTPUTS:
%   train_data:             (n x m x p x q) double; 4D array; 
%                           image x training_set_size x num_trials
%   train_labels:           (q x 1) double; Training label associated with
%                           each trial in train_data
%
% See also LOAD_DOPPLER_DATA

%% Input parser
p = inputParser;
p.addOptional('zscore',true, @islogical)
p.addOptional('spatial_filter',{'', [], []})
p.addOptional('training_set_size', 3);
p.addOptional('buffer_size', 60);
parse(p, varargin{:});
inputs = p.Results;


%% Extract variables from the data struct
behavior = data.behavior;
dop = data.dop;
timestamps = data.timestamps;

[yPix, xPix, ~] = size(dop);

%% Go through behavior, trial-by-trial and extract the relevant training data
% Interested in last frames of memory period.
target_acquire_start = parse_behavior_time(behavior, 'target_acquire');

valid_trials = find(~isnan(target_acquire_start));
num_valid_trials = length(valid_trials);
train_data = NaN(yPix, xPix, inputs.training_set_size, num_valid_trials);

train_labels = NaN(num_valid_trials, 1);

% Create progressbar
pb = waitbar(0, 'processing data');

% Iterate over each trial
for trial = 1:num_valid_trials
    % Find which frames we want
    trial_num = valid_trials(trial);
    current_trial_target_acquire_time = target_acquire_start(trial_num);
    last_frames_before_target_acquire = find(timestamps(:, 1) <= current_trial_target_acquire_time, inputs.training_set_size, 'last');
    
    % Create buffer
    begin_buff = max([1 last_frames_before_target_acquire(inputs.training_set_size)-inputs.buffer_size]);
    end_buff = last_frames_before_target_acquire(end);
    
    % Create z-scored data buffer
    data_buff = dop(:, :, begin_buff:end_buff);         
    data_buff = preprocess_data(data_buff,...
        'zscore',inputs.zscore,...
        'spatial_filter', inputs.spatial_filter);
    
    % Extract training data
    train_data(:, :, :, trial) = data_buff(:, :, end-inputs.training_set_size+1:end);
    
    % Assign label
    if any(data.actual_labels(last_frames_before_target_acquire) ~= data.actual_labels(last_frames_before_target_acquire(1)))
        error('Error in assigning training labels.');
    else
        train_labels(trial) = data.actual_labels(last_frames_before_target_acquire(1));
    end
    
    % Advance waitbar
    waitbar(trial/num_valid_trials);
end

% Close progressbar
close(pb);
end