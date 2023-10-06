function predict_movement_direction(data)
% predict_movement_direction  Essence of the real-time fUS-BMI
% This function takes in a data point and performs the fUS-BMI task using
% its information.
%
% Stripped version for simulated data for real-time fUS-BMI paper.
% Removed TCP socket information
%
%
% INPUTS:
%   data:           (struct); Multiple fields
%                   From real-time data:
%                       data - Power Doppler image
%                       ID   - Acquisition number
%                       time - Acquisition time
%                   For offline simulation, need the following additional
%                   fields
%                       class         - Indictes what direction or class
%                       phaseIn       - Indicates the state number
%                       nClasses      - Indicates how many classes.
%                       full_display  - Boolean; determine whether
%                                       to plot full display panel
%                       retrain       - Boolean to retrain the classifiers
%                                       (or just predict using the static
%                                       model)
%                       mask          - The boolean mask of which voxels to
%                                       use
%                       mem_length    - How many frames are in memory
%                                       period for this trial
%
% OUTPUTS:
%   None currently. Uses putvar instead.
%
% See also simulate_real_time_fUS_BMI


%% initializing persistents

% Specify the common persistents
persistent k Dop t_elapsed t label phase trial train...
    trainLabels class_predicted class_true model buffer_windows ...
    data_savepath num_frames_in_memory

% Specify the task_specific persistents
persistent class_predicted_horz class_predicted_vert class_true_horz class_true_vert ...
    custom_counting_matrix lookup_table_for_horz lookup_table_for_vert ...
    model_horz model_vert

%% Specifying task and parameters to run

% Approximate frame rate of the acquisition system
framerate = 2; % In Hz

% How many timepoints should we use from each trial. Relative to end of
% memory period
training_set_size = 3;

%% Passable arguments
full_display = data.full_display;
retrain = data.retrain;
use_previous_data_session = data.use_previous_data_session;
save_data = data.save_data;


nClasses = data.nClasses;
switch nClasses
    case 2
        task_type = '2 tgt';
    case 8
        task_type = '8 tgt';
    otherwise
        error('This number of classes is not supported.');
end

nTrials_before_train = data.nTrials_before_train;

filter_type = data.filter_type;
filter_size = data.filter_size;
filter_sigma = data.filter_sigma;

nb_block = data.num_frames;
class = data.class;
phase_in = data.phaseIn;
num_frames_in_memory = floor(data.mem_length);

if isempty(k)
    initialize_workspace()
end

%% update persistent values on each acquisition
k=k+1;
t(k) = now;
t_elapsed(k) = etime(datevec(now), datevec(t(1)));
Dop(:,:,k) = data.data;

if phase_in == 2
    % Reset for each trial
    num_frames_in_memory = [];
end

switch task_type
    case '8 tgt'
        if ~isnan(class)
            
            label(k) = class;
            
            label_horz = lookup_table_for_horz(class);
            label_vert = lookup_table_for_vert(class);
        else
            label(k) = NaN;
        end
    case '2 tgt'
        label(k) = class;
end
phase(k) = phase_in;

% if we're in a new phase & it's phase 2, increment trial counter
if phase(k)==2 && phase(k-1)~=2
    trial(k) = trial(k-1)+1;
else
    trial(k) = trial(k-1);
end

%%  separate the data in the buffer (using windows)
t_ind = false(size(t_elapsed));
start_ind = max(1, length(t_ind)-buffer_windows);
t_ind(start_ind:end) = true;

%% preprocess data (only used on this phase combination)
% Phase 4 = Memory
% Prepare for the testing set

if phase(k) == 4 && ...
        length(phase) > num_frames_in_memory &&...
        all(ismember(phase(k-num_frames_in_memory+1:k), 4)) &&...
        all(~ismember(phase(k-num_frames_in_memory), 4)) && ...
        size(train,1) > nTrials_before_train
    data_buff = Dop(:, :, t_ind);           % data in buffer
    data_buff = preprocess_data(data_buff,...
        'zscore',true,...
        'spatial_filter', {filter_type, filter_size, filter_sigma});
    [m, n, ~] = size(data_buff);
end

% Phase 8 = Success ITI. Prepare for the training set.
if length(phase) > 2 && ismember(phase(k), [8 9]) && ~ismember(phase(k-1), [8 9])
    data_buff = Dop(:, :, t_ind);           % data in buffer
    data_buff = preprocess_data(data_buff,...
        'zscore',true,...
        'spatial_filter', {filter_type, filter_size, filter_sigma});
    [m, n, ~] = size(data_buff);
end

% display the latest images coming into the function & buffer
if full_display || k==nb_block || nb_block==0
    plot_frame();
end


%% BCI
%%%%% ONLINE PREDICTION %%%%%
current_prediction = [];

% This defines the conditions for making a prediction.
if phase(k) == 4 && ...
        length(phase) > num_frames_in_memory &&...
        all(ismember(phase(k-num_frames_in_memory+1:k), 4)) &&...
        all(~ismember(phase(k-num_frames_in_memory), 4)) && ...
        size(train,1) > nTrials_before_train
    
    % create the test set
    new_test = reshape(data_buff(:,:,end-training_set_size+1:end), m, n*training_set_size);
    test = reshape(new_test, 1, []);
    
    % make the prediction using `train_classifier` and
    % `make_prediction`
    switch task_type
        case '8 tgt'
            % If no model has been trained yet, then train the model first.
            if isempty(model_horz) || isempty(model_vert)
                % Define the class labels for horizontal and vertical
                trainLabels_horz = lookup_table_for_horz(trainLabels);
                trainLabels_vert = lookup_table_for_vert(trainLabels);
                % Update the classifiers
                model_horz = train_classifier(train, ...
                    trainLabels_horz, ...
                    'method', 'PCA+LDA');
                model_vert = train_classifier(train, ...
                    trainLabels_vert, ...
                    'method', 'PCA+LDA');
            end
            % Make prediction
            
            class_horz = make_prediction(test, model_horz);
            class_vert = make_prediction(test, model_vert);
            
            
            % Combine horizontal and vertical prediction to make single
            % predicted class for 8-directions.
            current_prediction = class_horz + 3 * (class_vert-1);
            
            % Add prediction to log.
            class_predicted_horz = cat(2, class_predicted_horz, class_horz);
            class_predicted_vert = cat(2, class_predicted_vert, class_vert);
            class_predicted = cat(2, class_predicted, current_prediction);
            
            % Add prediction to confusion matrix
            custom_counting_matrix(class_predicted(end), label(k)) = custom_counting_matrix(class_predicted(end), label(k)) + 1;
            
            class_true_horz = cat(2, class_true_horz, label_horz);
            class_true_vert = cat(2, class_true_vert, label_vert);
            
        case '2 tgt'
            % If no model has been trained yet, then train the model first.
            if isempty(model) % update the classifiers & predict
                model = train_classifier(train, trainLabels, 'method', 'CPCA+LDA');
            end
            
            % Make prediction
            current_prediction = make_prediction(test, model);
            
            class_predicted= cat(2, class_predicted, current_prediction);
    end
    
    
    % save the actual class for later comparison
    class_true = cat(2, class_true, label(k));
    
    % image classification accuracy progress
    % plot_accuracy()   % plots accuracy as f(trial) vs. significance
    if full_display
        plot_counter();
        if nClasses == 8
            plot_split_performance();
            plot_confusion_matrix();
        end
    end
    
end % If statement for enough trials




% Plot performance in offline mode.
if (k==nb_block) || nb_block==0 || mod(k, 500) == 0 && size(train, 1) > nTrials_before_train && ~isempty(class_predicted)
    try
        plot_counter();
        if nClasses == 8
            plot_split_performance();
            plot_confusion_matrix();
        end
    catch
        warning('Unknown nonfatal error with plotting. Debug to fix.');
    end
end


%%%%% UPDATING TRAINING %%%%%
% if we have gotten through enough trials that classification might
% work, then let's just look at all previous trials' data as the
% training set & the current trial as the testing set.

if data.add_all_trials_to_training_set
    if length(phase) > 2 &&...
            ((ismember(phase(k), 8) && ~ismember(phase(k-1), 8)) || ...
            (ismember(phase(k), 9) && ismember(phase(k-1), [5, 6]))) &&...
            any(ismember(phase, [2, 3])) % This `and` statement handles possibility that we started in middle of trial on trial 1, and prevents it from breaking.
        
        % We only want to update the training if a prediction was made,
        % i.e. got to the target acquisition state. We might be missing a
        % few trials still if the prediction was made before we got to the
        % target acquisition state and the monkey made an error before
        % target acquisition.
        
        memory_start_ind = find(ismember(phase, [2, 3]), 1, 'last')+1;
        memory_start_ind_relative_to_buffer_end = k-memory_start_ind;
        % if we just started a go phase, update the training set and
        % re-train the classifiers
        
        % Use the last n timepoints before the movement period
        training_window_start = memory_start_ind_relative_to_buffer_end - num_frames_in_memory + training_set_size;
        training_window_end = memory_start_ind_relative_to_buffer_end - num_frames_in_memory+1;
        
        
        new_train = reshape(data_buff(:,:,end-training_window_start:end-training_window_end),m,n*3);
        
        train = cat(1, train, reshape(new_train, 1, []));
        
        trainLabels(size(train,1)) = label(k-1);
        
        % Train the classifier
        % Have this code in two places within the script because under full BMI
        % control, this code cannot be reached without a successful trial. The
        % successful trial requires a previously trained model.
        if size(train,1) > nTrials_before_train
            switch task_type
                case '8 tgt'
                    % Define the class labels for horizontal and vertical
                    trainLabels_horz = lookup_table_for_horz(trainLabels);
                    trainLabels_vert = lookup_table_for_vert(trainLabels);
                    
                    % Update the classifiers
                    if retrain || isempty(model_horz) || isempty(model_vert) ||...
                            (~retrain && size(train,1) == nTrials_before_train)
                        model_horz = train_classifier(train, ...
                            trainLabels_horz, ...
                            'method', 'PCA+LDA');
                        model_vert = train_classifier(train, ...
                            trainLabels_vert, ...
                            'method', 'PCA+LDA');
                    end
                    
                    
                case '2 tgt'
                    % If two classes, then predict and store prediction.
                    if retrain || isempty(model) ||...
                            (~retrain && size(train,1) == nTrials_before_train)% update the classifiers & predict
                        model = train_classifier(train, trainLabels, 'method', 'CPCA+LDA');
                    end
            end
        end
    end
else % Only add successful trials
    if length(phase) > 2 &&...
            ismember(phase(k), 8) &&...
            ~ismember(phase(k-1), 8) &&...
            any(ismember(phase, [2, 3])) % This `and` statement handles possibility that we started in middle of trial on trial 1, and prevents it from breaking.
        
        % We only want to update the training if it was a successful trial,
        % i.e. hits 'success_iti' state
        
        memory_start_ind = find(ismember(phase, [2, 3]), 1, 'last')+1;
        memory_start_ind_relative_to_buffer_end = k-memory_start_ind;
        % if we just started a go phase, update the training set and
        % re-train the classifiers
        
        % Use the last n timepoints before the movement period
        training_window_start = memory_start_ind_relative_to_buffer_end - num_frames_in_memory + training_set_size;
        training_window_end = memory_start_ind_relative_to_buffer_end - num_frames_in_memory+1;
        
        new_train = reshape(data_buff(:,:,end-training_window_start:end-training_window_end),m,n*training_set_size);
        
        train = cat(1, train, reshape(new_train, 1, []));
        
        trainLabels(size(train,1)) = label(k-1);
        
        % Train the classifier
        % Have this code in two places within the script because under full BMI
        % control, this code cannot be reached without a successful trial. The
        % successful trial requires a previously trained model.
        if size(train,1) > nTrials_before_train
            switch task_type
                case '8 tgt'
                    % Define the class labels for horizontal and vertical
                    trainLabels_horz = lookup_table_for_horz(trainLabels);
                    trainLabels_vert = lookup_table_for_vert(trainLabels);
                    
                    % Update the classifiers
                    if retrain || isempty(model_horz) || isempty(model_vert) ||...
                            (~retrain && size(train,1) == nTrials_before_train)
                        model_horz = train_classifier(train, ...
                            trainLabels_horz, ...
                            'method', 'PCA+LDA');
                        model_vert = train_classifier(train, ...
                            trainLabels_vert, ...
                            'method', 'PCA+LDA');
                    end
                    
                    
                case '2 tgt'
                    % If two classes, then predict and store prediction.
                    if retrain || isempty(model) ||...
                            (~retrain && size(train,1) == nTrials_before_train)% update the classifiers & predict
                        model = train_classifier(train, trainLabels, 'method', 'CPCA+LDA');
                    end
            end
        end
    end
end

%% data logger

if save_data
    if isempty(current_prediction)
        predicted_label = NaN;
    else
        predicted_label = current_prediction;
    end
    actual_label = label(k);
    phaseOut = phase(k);
    
    % Save some things in separate file for easier analysis later
    fidresult = fopen(fullfile(data_savepath, 'summaryResults.txt'), 'a');
    [fUS_comp_time, behavior_comp_time] = deal(datestr(data.time_acq, 'HH-MM-SS.FFF'));
    acq_timestamp = datestr(data.time_acq, 'HH:MM:SS.FFF');
    
    if ~isempty(num_frames_in_memory)
        num_frames_in_memory_forSave = num_frames_in_memory;
    else
        num_frames_in_memory_forSave = 0;
    end
    fprintf(fidresult, '%s %s %s %d %d %d %d\r\n', fUS_comp_time, behavior_comp_time, acq_timestamp, phaseOut, num_frames_in_memory_forSave, actual_label, predicted_label);
    fclose(fidresult);
end


%% update the task display
% Only run if in offlineMode

if full_display || k==nb_block || nb_block==0
    plot_task(phase, class, class_predicted, 2);
    
end
if k==nb_block
    putvar(class_true, class_predicted, use_previous_data_session, ...
        retrain, train, trainLabels);
end


%% Initialize workspace
    function initialize_workspace()
        k=1;                    % acquisition number (int)
        t(k) = now;             % absolute time (time)
        t_elapsed(k) = 0;       % elapsed time in sec (float)
        Dop(:,:,k) = data.data;      % data in time series format
        label(k) = 0;       % class label (int)
        buffer_windows = 60;    % windows
        phase(k) = 1;           % phase number (int)
        trial(k) = 1;           % trial number (int)
        
        custom_counting_matrix = zeros(9);
        
        data_savepath = data.savepath;
        
        % Turn off a warning for plotting
        warning('off', 'MATLAB:legend:IgnoringExtraEntries');
        
        
        if save_data
            % If directory already exists, clear the directory. Otherwise,
            % create folder
            if isfolder(data_savepath)
                delete_files = questdlg('This simulated Doppler data already exists. Do you want to overwrite this data file?',...
                            'delete old files?','yes','no','yes');
                if strcmp(delete_files, 'yes')
                    delete(fullfile(data_savepath, '*'));
                else
                    error('Check old files at "%s" and move them if desired before rerunning this script.', data_savepath);
                end
            else
                status = mkdir(data_savepath);
                if ~status
                    error('Save folder could not be created');
                end
            end
            
            
            add_all_trials_to_training_set = data.add_all_trials_to_training_set;
            % Save the experiment parameters
            save(fullfile(data_savepath, 'experiment_parameters.mat'), 'task_type', ...
                'retrain', 'use_previous_data_session', ...
                'buffer_windows', 'framerate', 'nTrials_before_train', ...
                'training_set_size', ...
                'add_all_trials_to_training_set', 'filter_type', ...
                'filter_size', 'filter_sigma');
        end
        
        % Initialize a dictionary
        % Create Container (equivalent to python dictionary)
        % These state names are the ones in Python. They won't perfectly
        % match with the state names from MATLAB.
        % These variables are not used, but keeping for clarity.
        keySet = {'init_fixation_acquire', 'init_fixation_hold', 'cue', 'memory', 'target_acquire', 'target_hold', 'reward', 'success_iti', 'abort_iti'};
        valueSet = [1                       2                      3       4           5               6               7       8              9];
        
        
        % create some lookup tables
        lookup_table_for_horz = [1 2 3 1 2 3 1 2 3];
        lookup_table_for_vert = [1 1 1 2 2 2 3 3 3];
        
        % This variable is not used, but keeping for clarity.
        lookup_table_8tgt = [5*pi/4 3*pi/2 7*pi/4 pi NaN 0 3*pi/4 pi/2 pi/4];
        
        
        if use_previous_data_session
            likely_loc = get_data_path('path_type', 'aligned');
            title_string = 'Find file used for the pre-training data';
            disp(title_string);
            [file, path] = uigetfile(fullfile(likely_loc, '*_training_data_full.mat'),...
                title_string);
            
            % If file, path is empty, then ask user if they want to create
            % an aligned dataset.
            if any([file path])
                prealigned_file = fullfile(path, file);
            else
                mode = questdlg('You did not select a file in the previous dialog box. Do you want to create an aligned dataset now?',...
                    'align new dataset?','yes','no','yes');
                if strcmp(mode, 'yes')
                    prealigned_file = create_aligned_training_dataset(...
                        'session_run_list', data.session_run_list);
                else
                    error('You specified you wanted to pretrain model but did not select or create an aligned dataset.');
                end          
            end
            
            
            load(prealigned_file, 'train_single', 'train_labels_single', ...
                'session_run_list_old');
            putvar(prealigned_file, session_run_list_old);
            fprintf('Previous training session loaded\n');
            train = double(train_single); % Convert back to double
            trainLabels = double(train_labels_single);
        end
        
        % Add some warnings
        if ~retrain
            warning('Retraining is not turned on. Are you sure this is what you wanted?');
        end
        if ~save_data
            warning('Data is not being saved. Make sure this is what you wanted!');
        end
        if use_previous_data_session
            warning('Using previous data session for training set. Make sure you have already run the alignment script');
        end
    end

%% plotting functions % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Because these are nested functions, they can access variables in the main
% function. Because of this, they cannot be easily moved into separate
% files...

    function plot_frame(subplot_start)
        
        if ~exist('subplot_start','var'), subplot_start = 1; end
        
        % create x & z inds
        z_mm = (0:size(data.data,1)-1)*0.1+3;
        x_mm = (0:size(data.data,2)-1)*0.1;
        
        % display the latest image
        set(gcf,'color','black')
        set(gcf, 'Units', 'normalized');
        set(gcf, 'Position', [0.2 0.2 0.7 0.7]);
        ax1 = subplot(2,3,subplot_start);
        colormap inferno;
        imagesc(x_mm,z_mm,nthroot(data.data,3));
        axis image
        colormap(ax1,inferno)
        %caxis(nthroot([0 1e+12],3))
        axis image
        
        % type set
        duration_string = duration(seconds(k/2),'format','hh:mm:ss.S');
        title(sprintf('trial time: %s',duration_string),'color','white')
        xticks(0:2:x_mm(end))
        yticks(z_mm(1):2:z_mm(end))
        ylabel('depth (mm)','color','white')
        xlabel('mm','color','white')
        set(ax1,'YColor','white','XColor','white')
        
        drawnow;
    end

    function plot_confusion_matrix(subplot_start)
        if ~exist('subplot_start','var')
            subplot_start = 6;
        end
        subplot(2, 3, subplot_start);
        set(gcf,'color','black');
        set(gcf, 'Units', 'normalized');
        set(gcf, 'Position', [0.2 0.2 0.7 0.7]);
        class_names = {'Left Down', 'Center Down', 'Right Down', ...
            'Left Center',  'Middle', 'Right Center', ...
            'Left Up', 'Center Up', 'Right Up'};
        
        class_reordering = [7 8 9 6 3 2 1 4 5];
        reordered_label_strings = class_names(class_reordering);
        
        reordered_counting_matrix = custom_counting_matrix(class_reordering, class_reordering)';
        
        % #Add this function or change
        confusion_chart = create_confusion_matrix(reordered_counting_matrix(1:8, :), ...
            'label_strings', reordered_label_strings, ...
            'title', 'Combined prediction', ...
            'FontColor', 'w');
        
        axis image;
        
        
        
        % Plot confusion matrix
        drawnow;
    end

    function plot_split_performance(subplot_start)
        if ~exist('subplot_start', 'var')
            subplot_start = 4;
        end
        set(gcf,'color','black');
        set(gcf, 'Units', 'normalized');
        set(gcf, 'Position', [0.2 0.2 0.7 0.7]);
        n_trials = length(class_true);
        
        % get accuracy for horizontal
        correct_count_horz = sum(class_predicted_horz==class_true_horz);
        correct_percent_horz = correct_count_horz/n_trials*100;
        
        % get accuracy for vertical
        correct_count_vert = sum(class_predicted_vert==class_true_vert);
        correct_percent_vert = correct_count_vert/n_trials*100;
        
        % Find threshold for binomial test
        multicoder_classes = 3;
        sig_thresh = binoinv(0.95, n_trials, 1/multicoder_classes);
        
        % make plot for horizontal
        ax4 = subplot(2,3,subplot_start);
        cla
        plot([0 3], [sig_thresh sig_thresh],'--r','linewidth',2)
        plot([0 3], [n_trials/multicoder_classes n_trials/multicoder_classes],'--k','linewidth',1)
        bar(correct_count_horz,'k'), hold on
        axis([0 2 0 n_trials*1.5])
        
        % type set
        legend('binomial test, p=0.05','chance',...
            'location','northwest')
        legend('boxoff')
        yticks('auto')
        new_yticks = unique(round(yticks));
        new_yticks(new_yticks>n_trials) = [];
        yticks(new_yticks)
        ylabel('number of trials')
        xticks(1);
        xticklabels({'correct'});
        set(ax4,'xcolor',[1 1 1]), set(ax4,'ycolor',[1 1 1])
        
        % print correct/incorrect percentages on figure
        text(1, correct_count_horz*1.01, sprintf('%2.1f%%',correct_percent_horz),...
            'VerticalAlignment','bottom','HorizontalAlignment','center')
        title('Horizontal accuracy', ...
            'Color', 'w');
        
        % make plot for vertical
        ax5 = subplot(2,3,subplot_start+1);
        cla
        plot([0 3], [sig_thresh sig_thresh],'--r','linewidth',2)
        plot([0 3], [n_trials/multicoder_classes n_trials/multicoder_classes],'--k','linewidth',1)
        bar(correct_count_vert,'k'), hold on
        axis([0 2 0 n_trials*1.5])
        
        % type set
        legend('binomial test, p=0.05','chance',...
            'location','northwest')
        legend('boxoff')
        yticks('auto')
        new_yticks = unique(round(yticks));
        new_yticks(new_yticks>n_trials) = [];
        yticks(new_yticks)
        ylabel('number of trials')
        xticks(1);
        xticklabels({'correct'});
        set(ax5,'xcolor',[1 1 1]), set(ax5,'ycolor',[1 1 1])
        
        % print correct/incorrect percentages on figure
        text(1, correct_count_vert*1.01, sprintf('%2.1f%%',correct_percent_vert),...
            'VerticalAlignment','bottom','HorizontalAlignment','center')
        title('Vertical accuracy', ...
            'Color', 'w');
        
        
        
        drawnow;
    end

    function plot_counter(subplot_start)
        
        if ~exist('subplot_start','var'), subplot_start = 3; end
        set(gcf,'color','black');
        set(gcf, 'Units', 'normalized');
        set(gcf, 'Position', [0.2 0.2 0.7 0.7]);
        % print the result
        n_trials = length(class_true);
        
        % get accuracy
        correct_count = sum(class_predicted==class_true);
        correct_percent = correct_count/n_trials*100;
        
        fprintf('%i/%i predicted correctly (%3.2f %%)\n', ...
            correct_count, n_trials, correct_percent);
        
        % get confidence intervals
        sig_thresh = binoinv(0.95, n_trials, 1/data.nClasses)+1;
        
        % make plot
        ax3 = subplot(2,3,subplot_start);
        cla
        plot([0 3], [sig_thresh sig_thresh],'--r','linewidth',2)
        plot([0 3], [n_trials/data.nClasses n_trials/data.nClasses],'--k','linewidth',1)
        bar(correct_count,'k'), hold on
        axis([0 2 0 n_trials*1.5])
        
        % type set
        legend('binomial test, p=0.05','chance',...
            'location','northwest')
        legend('boxoff')
        yticks('auto')
        new_yticks = unique(round(yticks));
        new_yticks(new_yticks>n_trials) = [];
        yticks(new_yticks)
        ylabel('number of trials')
        xticks(1);
        xticklabels({'correct'});
        set(ax3,'xcolor',[1 1 1]), set(ax3,'ycolor',[1 1 1])
        
        % print correct/incorrect percentages on figure
        text(1, correct_count*1.01, sprintf('%2.1f%%',correct_percent),...
            'VerticalAlignment','bottom','HorizontalAlignment','center')
        title('Overall accuracy', ...
            'Color', 'w');
        
        drawnow;
    end

    function plot_task(phase, class, predicted, subplot_start)
        
        % TargetPos for each class
        targetPos_8 = [19 0;
            13.43 13.43;
            0 19;
            -13.43 13.43;
            -19 0;
            -13.43 -13.43;
            0 -19
            13.43 -13.43];
        
        targetPos_9 = [-13.43 -13.43;
            0 -19;
            13.43 -13.43;
            -19 0;
            0 0;
            19 0;
            -13.43 13.43
            0 19
            13.43 13.43];
        
        targetPos_2 = [-19 0;
            19 0;];
        
        predicted_lookup_table = [targetPos_8; 0 0;];
        predicted_lookup_table = predicted_lookup_table([6 7 8 5 9 1 4 3 2], :);
        
        % set up the axis
        if ~exist('subplot_start','var'), subplot_start = 2; end
        ax2 = subplot(2,3,subplot_start);
        set(gcf,'color','black');
        set(gcf, 'Units', 'normalized');
        set(gcf, 'Position', [0.2 0.2 0.7 0.7]);
        if ~isempty(predicted)
            predicted = predicted(end);
        end
        
        switch phase(end)
            case 0 %Trial start - Not reliably acquired
                create_task_screen(NaN, NaN, NaN, 'trial start', ...
                    'show_center_fixation', false);
            case 1 %Fixation acquire - Not reliably acquired
                create_task_screen(NaN, NaN, NaN, 'fixation acq');
            case 2 %Fixation hold
                create_task_screen(NaN, NaN, NaN, 'fixation hold');
            case 3 %Cue - Not reliably acquired
                switch nClasses
                    case 2
                        create_task_screen(targetPos_2(class, :), NaN, NaN, 'cue');
                    case 8
                        create_task_screen(targetPos_9(class, :), NaN, NaN, 'cue');
                end
            case 4  % memory
                % first frame in memory phase is a directional cue
                if phase(end-1)~=4
                    switch nClasses
                        case 2
                            create_task_screen(targetPos_2(class, :), NaN, NaN, 'cue');
                        case 8
                            create_task_screen(targetPos_9(class, :), NaN, NaN, 'cue');
                    end
                else
                    create_task_screen(NaN, NaN, NaN, 'memory');
                end
            case 5 % Target acquire
                % Not reliably acquired
            case 6 % Target hold
                % Not reliably acquired (but more reliable than some other
                % phases)
            case 7 % Reward
                % Not reliably acquired
                
            case {8, 9}  % iti
                if ~ismember(phase(end-1), [8, 9]) % Use first time point as move
                    % training
                    switch nClasses
                        case 2
                            if isempty(predicted)
                                create_task_screen(targetPos_2(class, :), targetPos_2(class, :), NaN, 'move', ...
                                    'show_center_fixation', false);
                                % testing
                            else
                                create_task_screen(targetPos_2(class, :), NaN, targetPos_2(predicted, :), 'predict', ...
                                    'show_center_fixation', false);
                            end
                        case 8
                            if isempty(predicted)
                                create_task_screen(targetPos_9(class, :), targetPos_9(class, :), NaN, 'move', ...
                                    'show_center_fixation', false);
                            else
                                create_task_screen(targetPos_9(class, :), NaN, predicted_lookup_table(predicted, :), 'predict', ...
                                    'show_center_fixation', false);
                            end
                    end
                elseif ~ismember(phase(end-2), [8, 9]) && ismember(phase(end-1), [8 9])
                    switch nClasses
                        case 2
                            % training period
                            if isempty(predicted)
                                create_task_screen(targetPos_2(class, :), targetPos_2(class, :), NaN, 'hold', ...
                                    'show_center_fixation', false);
                                % testing period
                            else
                                create_task_screen(targetPos_2(class, :), targetPos_2(class, :), targetPos_2(predicted, :), 'hold', ...
                                    'show_center_fixation', false);
                            end
                        case 8
                            if isempty(predicted)
                                create_task_screen(targetPos_9(class, :), targetPos_9(class, :), NaN, 'hold', ...
                                    'show_center_fixation', false);
                            else
                                create_task_screen(targetPos_9(class, :), targetPos_9(class, :), predicted_lookup_table(predicted, :), 'hold', ...
                                    'show_center_fixation', false);
                            end
                    end
                else
                    create_task_screen(NaN, NaN, NaN, 'iti', ...
                        'show_center_fixation', false);
                end
            otherwise
                %fprintf('Phase %d \n', phase(end));
                
        end
        xlabel(sprintf('trial #%i',trial(end)),'color','white')
        
        
        drawnow;
    end

end
