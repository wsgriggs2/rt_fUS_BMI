%% Script that will simulate real-time fUS-BMI
% Loads in data, simulates streamining data in to the decoder frame by
% frame, and performs decoding.

%% Initialize workspace
clear; close all; clc;

%% Load project record for reference
% Need to include the simulated results here so that we can properly
% increment the run count.
project_record = get_project_record('include_simulated', true);


%% select run mode
fast_mode = 'Run as fast as possible with minimal user display';
verbose_mode = 'Run slower with full display of data streaming and trial performance';

mode = questdlg('Which mode would you like to run in?',...
    'Mode selection',...
    fast_mode, ...
    verbose_mode,...
    fast_mode);

% Dialog box to select desired settings
option_strings = {'pretrain+retrain', 'pretrain only', 'retrain only'};
retrain_and_pretrain_choice = listdlg('ListString', option_strings, ...
    'PromptString','Select desired combination of training methods:', ...
    'SelectionMode', 'single');
switch option_strings{retrain_and_pretrain_choice}
    case 'pretrain+retrain'
        retrain = true;
        use_previous_data_session = true;
    case 'retrain only'
        retrain = true;
        use_previous_data_session = false;
    case 'pretrain only'
        retrain = false;
        use_previous_data_session = true;
end

% Some standard parameters
nTrials_before_train = 40; 
add_all_trials_to_trainingset = true;
save_data = true;

% Define spatial filter parameters
filter_to_use = 'disk';
filter_size = 2;
filter_sigma = 0;

%% load data
data = load_doppler_data;

% Check to make sure only single session is loaded
if size(data.session_run_list, 1) > 1
    error('This script is only designed to work with a single session being loaded.');
end

behavior = data.behavior;
timestamps = data.timestamps;
dop = data.dop;
coreParams = data.coreParams;

%Specify the save path if we are saving data
if save_data
    session_num = data.session_run_list(:,1);
    run_num = data.session_run_list(:,2);
    indx = find(ismember(project_record.Session, session_num) & ismember(project_record.Run, run_num));

    % Increment run number
    session_indx = ismember(project_record.Session, project_record.Session(indx));
    run_numbers = project_record.Run(session_indx);
    next_run_number = max(run_numbers) + 1;
    
    
    data_savepath = fullfile(get_data_path('path_type', 'simulated'), ...
        sprintf('S%d_R%d', session_num, next_run_number));
else
    data_savepath = '';
end

%% Parse Behavior
% Assigning labels to multi-targets
targetPos = vertcat(behavior.targetPos);
if iscell(targetPos)
    targetPos = cell2mat(targetPos);
end
[angle,distance] = cart2pol(targetPos(:,1),targetPos(:,2)); % convert to polar coordinates
angle(angle<0) = angle(angle<0)+2*pi; %Convert to have only positive angles.

UniqueAngles = unique(angle);
TargetPosInd = zeros(length(targetPos),1);
for position = 1:length(UniqueAngles)
    TargetPosInd(angle==UniqueAngles(position)) = position;
end

switch length(UniqueAngles)
    case 8
        decode_type = '8 tgt'; % Options are '2 tgt' and '8 tgt';
        nClasses = 8;
    case 2
        decode_type = '2 tgt'; % Options are '2 tgt' and '8 tgt';
        nClasses = 2;
end

switch decode_type
    case '8 tgt'
        train_labels = NaN(size(targetPos,1),2);
        
        figure;
        tiledlayout(1,2);
        
        for dimension = 1:2
            train_labels(targetPos(:, dimension) < 0, dimension) = 1;
            
            train_labels(targetPos(:, dimension) == 0, dimension) = 2;
            train_labels(targetPos(:, dimension) > 0, dimension) = 3;
            
            nexttile;
            histogram(train_labels(:,dimension));
            title(sprintf('Number of trials for each training label: dimension %d',dimension));
        end
        
        behavior_trim = behavior;
        
    case '2 tgt'
        if isequal(UniqueAngles, [0:pi/4:2*pi-pi/16]')
            % If this is from an 8 target dataset, then we don't want the
            % up/down conditions. We will keep the targets that are
            % contra/ipsi though. Can also restrict to just a subset, such
            % as 1, 5 (pure L/R), or 2, 4 (upwards L/R), or 6, 8 (downwards
            % L/R).
            behavior_trim = behavior;
            discard_ind = ismember(TargetPosInd, [3, 7]);
            behavior_trim(discard_ind) = [];
        elseif isequal(UniqueAngles, [0; pi])
            behavior_trim = behavior;
        else
            error('Cannot handle this dataset. Unknown error');
        end
    otherwise
        error('This decode_type is not supported');
end

%% get important data & task variables
[~, ~, n_frames] = size(dop);   % get data size

% rebuilding task structure according to real-time BMI convention:
task_labels = {'trialstart', 'initialfixation', 'fixationhold', 'cue', ...
    'memory', 'target_acquire', 'target_hold', 'reward', 'iti'};

%% set up video save

% ask the user if they want to save
if strcmp(mode, verbose_mode)
    save_resp = questdlg('Would you like to save the data & video?',...
        'save','yes','no','yes');
    save_bool = strcmp(save_resp,'yes');
else
    save_bool = false;
end

% if they do, set up the video object for writing
if save_bool
    title_string = 'Specify where to save video';
    disp(title_string);
    video_filepath = uigetdir(get_data_path('path_type', 'output'), title_string);
    video_filename = strcat('realtime_BCI_panel_',...
        datestr(date,'yyyy-mm-dd'));
    video_filename = fullfile(video_filepath, video_filename);
    vidObj = VideoWriter(video_filename,'MPEG-4');
    vidObj.FrameRate = 24;
    open(vidObj);
end

%% run simulation
clear predict_movement_direction

% set up simulation timing
tic();
t_elapsed = toc();
train_time = [];
predict_time = [];
trial = 1;

[phase_num, class] = deal(zeros(n_frames, 1));

% set up the figure
close all;
f1 = figure(1);
if strcmp(mode, verbose_mode)
    set(f1, 'Units', 'normalized', 'position', [0.1 0.1 0.9 0.9])
end
for frame = 1:n_frames
    % Build struct to be passed to predict_movement_direction
    data_struct.data = dop(:, :, frame);
    data_struct.time_acq = data.time_acq(frame);
    data_struct.phaseIn = data.state_num(frame);
    data_struct.class = data.actual_labels(frame);
    data_struct.ID = frame;
    data_struct.time = timestamps(frame, 1);
    data_struct.pos = [];
    data_struct.full_display = false;
    data_struct.retrain = retrain;
    data_struct.mem_length = data.memory_length(frame);
    data_struct.num_frames = n_frames;
    data_struct.use_previous_data_session = use_previous_data_session;
    data_struct.nTrials_before_train = nTrials_before_train;
    data_struct.add_all_trials_to_training_set = add_all_trials_to_trainingset;
    data_struct.save_data = save_data;
    data_struct.savepath = data_savepath;
    data_struct.filter_type = filter_to_use;
    data_struct.filter_size = filter_size;
    data_struct.filter_sigma = filter_sigma;
    data_struct.session_run_list = data.session_run_list;
    
    switch decode_type
        case '8 tgt'
            data_struct.nClasses = 8;
        case '2 tgt'
            data_struct.nClasses = 2;
    end
    
    % print time & phase info to screen to screen
    if strcmp(mode, verbose_mode)
        % full display mode calls predict_movement_direction & prints full display panel
        data_struct.full_display = true;
        predict_movement_direction(data_struct)
    elseif strcmp(mode, fast_mode)
        % time mode runs predict_movement_direction as fast as possible (no panels)
        tmp = tic();
        predict_movement_direction(data_struct)
        predict_time = cat(2, predict_time, toc(tmp));
    end
    
    % update the timing
    t_elapsed(frame) = toc();
    
    % update the trial & frame we'll use on the next loop
    if frame > 1 && data.state_num(frame-1) > data.state_num(frame)
        trial = trial+1;
    end
    
    % update the frame & write the frame to video
    if strcmp(mode, verbose_mode)
        drawnow
    end
    if strcmp(mode, verbose_mode) && save_bool
        drawnow
        writeVideo(vidObj,getframe(gcf));
    end
end

% close the video object
if strcmp(mode, verbose_mode)
    close(vidObj)
end

%% Update project record
if save_data
    
    new_entry = project_record(indx, :);
    
    % Find maximum run number for the session. Start incrementally using
    % run numbers.
    
    % Need to write in way to load both session records and combine them to
    % get an accurate run #.
    session_indx = ismember(project_record.Session, project_record.Session(indx));
    run_numbers = project_record.Run(session_indx);
    new_entry.Run = next_run_number;
    new_entry.simulated_results = true;
    
    new_entry.pretraining = double(use_previous_data_session);
    new_entry.retraining = double(retrain);
    if new_entry.pretraining
        new_entry.pretraining_sessionrun = sprintf('S%dR%d', session_run_list_old(1, 1), session_run_list_old(1, 2));
        new_entry.pretraining_file = prealigned_file;
    else
        new_entry.pretraining_sessionrun = 'N/A';
        new_entry.pretraining_file = 'N/A';
    end
    
    new_entry.Notes = {sprintf('Simulated results for S%dR%d generated by "simulate_real_time_fUS_BMI.m" on %s.', data.session_run_list(1, 1), data.session_run_list(1, 2), datestr(now(), 'YYYY-mm-dd'))};
    
    % Add the new entry
    add_entry_to_project_record(new_entry);
    
end

%% Plot performance
if size(data.session_run_list, 1)>1
    session_string = sprintf('%d ', data.session_run_list(:, 1)');
    subtitle = sprintf('[%s ]\n Retrain: %s PreviousData: %s', session_string, string(retrain), string(use_previous_data_session));
else
    subtitle = sprintf('S%dR%d\n Retrain: %s PreviousData: %s', data.session, data.run, string(retrain), string(use_previous_data_session));
end

plot_performance(class_true, class_predicted, ...
    'subtitle', subtitle, ...
    'pvalue_threshold', 0.05, ...
    'n_training_trials', nTrials_before_train, ...
    'num_replicates', 1000, ...
    'rolling_average', 100);


