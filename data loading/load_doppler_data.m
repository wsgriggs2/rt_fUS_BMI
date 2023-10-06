function data_out = load_doppler_data(varargin)
% LOAD_DOPPLER_DATA  Load specified dataset
%
% INPUTS:
%   varargin: 
%       gui_title:              string; Title for GUI
%       session_run_list:       (1 x 2) double; [session run]' If specified, 
%                               skips selection GUI. Should match session/run
%                               numbers from `project_record.json`.
%       verbose:                bool; Do you want list of sessions/runs
%                               printed to command line?
%       data_path:              string; Path to data
%       manual_alignment:       bool; Align all specified sessions?
%       single_or_multiple:     string; {'single', 'multiple'}; Load one or
%                               multiple sessions?
%       variables_to_load:      (1 x n) cell array of strings; Names of the
%                               variables we want to load. Takes advantage
%                               of HDF5 data structure of .mat files to
%                               speed up loading if we do not need all
%                               variables.
%      
% OUTPUTS:
%   data_out:                   Struct; The specified loaded data. If
%                               multiple sessions were specified, they are 
%                               all concatenated rather than kept separate.



%% Input parser

possible_variables = {'dop', 'behavior', 'UF', 'timestamps', ...
    'neurovascular_map', 'actual_labels', 'predicted_labels', 'coreParams', ...
    'frame_num', 'memory_length', 'state_num', 'time_acq', ...
    'session', 'run'};

p = inputParser;
addOptional(p, 'gui_title', 'Select Session/Run(s) to load:');
addOptional(p, 'session_run_list', []);
addOptional(p, 'verbose', true, @islogical);
addOptional(p, 'data_path', get_data_path('path_type', 'doppler'));
addOptional(p, 'manual_alignment', true, @islogical);
addOptional(p, 'single_or_multiple', 'single'); % Can be either 'single' or 'multiple'
addOptional(p, 'variables_to_load', possible_variables);
parse(p, varargin{:});
inputs = p.Results;

% load Session Info
project_record = get_project_record('include_simulated', false);


%% Declaring some other needed variables
state_names = {'trialstart', 'initialfixation', 'fixationhold', 'cue', ...
    'memory', 'target_acquire', 'target_hold', 'reward', 'iti'};


%% Select the desired session
if isempty(inputs.session_run_list)
    % create the list of available sessions/runs
    runStrings = cell(size(project_record,1),1);
    
    for i = 1:size(project_record,1)
        % Extract information
        Session = project_record.Session(i);
        Run = project_record.Run(i);
        Slot = project_record.Slot(i);
        Monkey = project_record.Monkey(i);
        nTargets = project_record.nTargets(i);
        pretraining = logical(project_record.pretraining(i));
        retraining = logical(project_record.retraining(i));
        
        % Create string
        runStrings{i} = ['Monkey ' char(Monkey) ', Session ' num2str(Session) ...
            ', Run ' num2str(Run) ', Slot ' num2str(Slot) ...
            ', nTargets ' num2str(nTargets) ', pretrain:' mat2str(pretraining) ...
            ', retrain:' mat2str(retraining)];
        
    end
    
    % use a GUI to select which ones you want
    [indx,tf] = listdlg('PromptString', inputs.gui_title,...
        'SelectionMode',inputs.single_or_multiple,...
        'ListString',runStrings, ...
        'ListSize', [600 300]);
    
    
    session_num = project_record.Session(indx);
    run_num = project_record.Run(indx);
    session_run_list = [session_num run_num];
else
    session_run_list = inputs.session_run_list;
    session_num = session_run_list(:,1);
    run_num = session_run_list(:,2);
    indx = ismember(project_record.Session, session_num) & ismember(project_record.Run, run_num);
    tf = true;
end

if ~tf
    return;
end

%% Make sure specified sessions are all compatible

% Check if more than one monkey is used
monkey_pass = length(unique(project_record.Monkey(indx)))==1;

% Check if all are same effector
effector_pass = length(unique(project_record.Task(indx)))==1;

% Check if all are same num of targets
task_pass = length(unique(project_record.nTargets(indx)))==1;

% Check if all from same imaging plane
plane_pass = length(unique(project_record.Slot(indx))) == 1;

% Pause if the sessions are deemed incompatible
if ~(monkey_pass && task_pass && effector_pass && plane_pass)
    error('The selected sessions are not compatible and cannot be combined!');
end

%% Load the selected session
n_sessions_to_load = size(session_run_list, 1);
for session = 1:n_sessions_to_load
    
    current_session = session_num(session);
    current_run = run_num(session);
    
    % Define filename and filepath
    doppler_filename = sprintf('rt_fUS_data_S%d_R%d.mat', current_session, current_run);
    data_fullpath = fullfile(inputs.data_path, doppler_filename);
    
    % Load file
    data_in = load(data_fullpath);
    
    % Adjust timing in behavior
    [data_in.behavior, data_in.timestamps] = clean_behavior_and_timestamps(data_in.behavior, data_in.timestamps);
    
    % Add field to struct to represent session number and run number
    [data_in.behavior.session] = deal(current_session);
    [data_in.behavior.run] = deal(current_run);
     
    % Add fields to timestamps to represent session number and run number
    data_in.timestamps(:, 2) = current_session;
    data_in.timestamps(:, 3) = current_run;
    
    if session == 1
        % Initialize the output data array
        data_out = data_in;
    else
        %% Adjust timing so taht the timestamps do not overlap between
        % different sessions
        
        % Find the last timestamp in `data_out.behavior`
        last_behavior_entry = data_out.behavior(end);
        for k = length(state_names):-1:1
            possible_last_timestamp = last_behavior_entry.(state_names{k});
            if ~isempty(possible_last_timestamp)
                last_state_in_final_trial = state_names{k};
                break;
            end
        end
        
        % Apply timestamp correction
        maxTime = max([parse_behavior_time(data_out.behavior(end), last_state_in_final_trial), data_out.timestamps(:, 1)'], [], 'all');
        data_in.timestamps(:, 1) = data_in.timestamps(:, 1) + maxTime + 1; % Shift by 1 additional second to avoid overlap of any timestamps
        data_in.behavior = adjust_behavior_time(data_in.behavior, state_names, maxTime);
        
        %% Align session's anatomy to previous session(s)
        if inputs.manual_alignment
            
            anatomy_prev = data_out.neurovascular_map;
            anatomy_curr = data_in.neurovascular_map;
            
            % Use manual alignment to register between sessions/runs
            alignment_tform = align_neurovascular_maps(anatomy_prev, anatomy_curr);
            %Convert transform into matlab format
            % Alignment_tform is passed back from ManualImageAlignment
            tform = affine2d(alignment_tform);
            % Align image to match previously loaded data
            data_in.dop = imwarp(...
                data_in.dop, ...
                tform, ...
                'OutputView', imref2d(size(anatomy_prev)));
        end
        
        %% Make sure the data is the same size before combining
        if ~all([size(data_out.dop,1)==size(data_in.dop,1),size(data_out.dop,2)==size(data_in.dop,2)])
            % Use minimum spatial dimensions and minimum time dimension
            minDims = min(size(data_out.dop), size(data_in.dop));
            data_out.dop = data_out.dop(1:minDims(1), 1:minDims(2), :);
            data_in.dop = data_in.dop(1:minDims(1), 1:minDims(2), :);
        end
        
        % Concatenate dop in time
        data_out.dop =  cat(3, data_out.dop, data_in.dop);
        
        
        %% update UF
        if ~isfield(data_in.UF,'Depth')
            warning('UF missing depth field')
        elseif ~isequal(data_in.UF.Depth, data_out.UF.Depth)
            warning('Different depths between the sessions. Using smaller depth');
        end
        
        %% Concatenate behavior
        % Handle structs with different number of fields
        behaviorOut_fields = fieldnames(data_out.behavior);
        behavior_fields = fieldnames(data_in.behavior);
        if length(behaviorOut_fields) >  length(behavior_fields)
            warning('data_out.behavior has more fields than BehaviorIn. Make sure the new fields added make sense.');
            new_fields = behaviorOut_fields(~ismember(behaviorOut_fields, behavior_fields));
            empty_cells = cell(length(data_in.behavior), 1);
            for field = 1:length(new_fields)
                [data_in.behavior.(new_fields{field})] = empty_cells{:};
            end
        elseif length(behaviorOut_fields) <  length(behavior_fields)
            warning('data_out.behavior has more fields than BehaviorIn. Make sure the new fields added make sense.');
            behaviorOut_fields = fieldnames(data_out.behavior);
            behavior_fields = fieldnames(data_in.behavior);
            new_fields = behavior_fields(~ismember(behavior_fields, behaviorOut_fields));
            empty_cells = cell(length(data_out.behavior), 1);
            for field = 1:length(new_fields)
                [data_out.behavior.(new_fields{field})] = empty_cells{:};
            end
        end
        data_out.behavior = [data_out.behavior data_in.behavior];
        fprintf('Successful trials concatenated from %i to %i\n', ...
            size(data_out.behavior,2)-size(data_in.behavior,2)+1, size(data_out.behavior,2))
        
        
        %% Concatenate the timestamps
        data_out.timestamps = [data_out.timestamps; data_in.timestamps];
        
        %% Concatenate the BCI performance variables
        data_out.actual_labels = [data_out.actual_labels; data_in.actual_labels];
        data_out.predicted_labels = [data_out.predicted_labels; data_in.predicted_labels];
        data_out.state_num = [data_out.state_num; data_in.state_num];
        data_out.time_acq = [data_out.time_acq; data_in.time_acq];
        data_out.memory_length = [data_out.memory_length; data_in.memory_length];
        
        % Offset the new framenumbers
        data_out.frame_num = [data_out.frame_num; data_in.frame_num + data_out.frame_num(end)];  
    end % End if statement 
end % End for loop over each session

%% Some last cleaning of data loaded
data_out.session_run_list = session_run_list;
inputs.variables_to_load = [inputs.variables_to_load 'session_run_list'];

% This implementation does not take advantage of HDF5 structure of .mat
% file. Messy to implement due to concatenation but can specify which
% variables we want to load at the `load` function instead of loading all 
% variables, then removing the non-requested extraneous variables here.
possible_variables = fieldnames(data_out);
fields_to_remove = setdiff(possible_variables, inputs.variables_to_load);
if ~isempty(fields_to_remove)
    data_out = rmfield(data_out, fields_to_remove);
end
