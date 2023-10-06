%% Script to visualize data
% Works for both actual and simulated results


%% User settings
num_replicates = 1000;

update_confusion_matrix_for_each_trial = false;
pause_length = 0.05;

pvalue_threshold = 0.05;

rolling_window_size = NaN; % If NaN, then won't use a rolling widnow.


%% Define some useful things
combined_label_strings = {'Left Down', 'Center Down', 'Right Down', ...
    'Left Center',  'Middle', 'Right Center', ...
    'Left Up', 'Center Up', 'Right Up'};
reordered_label_strings = combined_label_strings([7 8 9 6 3 2 1 4 5]);


%% Specify files to be loaded
% load Session Info
project_record = get_project_record('include_simulated', true, ...
    'verbose', true);

% create the list of available sessions/runs
runStrings = cell(size(project_record,1),1);
for i = 1:size(project_record,1)
    Session = project_record.Session(i);
    Run = project_record.Run(i);
    Slot = project_record.Slot(i);
    Monkey = project_record.Monkey(i);
    nTargets = project_record.nTargets(i);
    retraining = project_record.retraining(i);
    pretraining = project_record.pretraining(i);
    runStrings{i} = ['Monkey ' char(Monkey) ', Session ' num2str(Session) ', Run ' num2str(Run) ', Slot ' num2str(Slot) ...
        ', nTargets ' num2str(nTargets)...
        ', retraining ' mat2str(retraining) ...
        ', prev session ' mat2str(pretraining)];
end

% use a GUI to select which ones you want
[indx,tf] = listdlg('PromptString','Select Session/Run(s) to load:',...
    'SelectionMode','multiple',...
    'ListString',runStrings, ...
    'ListSize', [600 300]);

session_run_list = [project_record.Session(indx) project_record.Run(indx)];

%% Define some useful variables
nTargets = project_record.nTargets(indx);
effector = project_record.Task(indx);
monkey = project_record.Monkey(indx);
slot = project_record.Slot(indx);
if any(nTargets ~= nTargets(1)) || any(monkey ~= monkey(1)) || ...
        any(slot ~= slot(1)) || any(effector ~= effector(1))
    error('The selected session are not compatible. Please rechoose your sessions to be analyzed.');
end

% Since they are all equal, we only need the first entry.
nTargets = nTargets(1);
monkey = monkey(1);
slot = slot(1);
effector = effector(1);

% How many sessions are we analyzing
num_sessions = size(session_run_list, 1);


%% Set up some things depending on nTargets
% Initialize variables
nTrials_for_plotting = []; 
accuracy_range_for_plotting = [];
angular_error_for_plotting = [];
accuracy_colorbar_limits = [];

% If analyzing more than one session, let's standardize appearances.
if num_sessions > 1
switch nTargets
    case 2
        switch effector
            case 'saccade'
                nTrials_for_plotting = [1 425]; %For 2-tgt saccade
                accuracy_range_for_plotting = [0 1];
                angular_error_for_plotting = [0 pi];
                accuracy_colorbar_limits = [];
            case 'reach'
                nTrials_for_plotting = [1 145]; %For 2-tgt reach
                accuracy_range_for_plotting = [0 1];
                angular_error_for_plotting = [0 pi];
                accuracy_colorbar_limits = [];
        end
    case 8
        nTrials_for_plotting = [1 362]; %For 8-tgt saccade
        accuracy_range_for_plotting = [0 1];
        angular_error_for_plotting = [0 3*pi/4];
        accuracy_colorbar_limits = [0 100];
end
end

%% Preallocate variables
switch nTargets
    case 2
        confusion_mat_acrossSession = NaN(2, 2, num_sessions);
    case 8
        confusion_mat_acrossSession = NaN(9, 9, num_sessions);
end
[empiric_mean_angular_error, null_mean_angular_error, pvalue_empiric, ...
    performance_acrossSession, predicted_class_noZeros, ...
    actual_class_noZeros] = deal(cell(num_sessions, 1));
[number_training_trials, last_nonsignificant_trial] = deal(NaN(num_sessions, 1));

%% Load necessary data
for session_counter = 1:num_sessions
    
    session = project_record.Session(indx(session_counter));
    run = project_record.Run(indx(session_counter));
    
    if project_record{indx(session_counter), 'simulated_results'}
        % Specify location of summary file for session
        base_path = get_data_path('path_type', 'simulated');
        subfolder_path = sprintf('S%d_R%d', session, run);
        summary_results_fname = fullfile(base_path, subfolder_path, 'summaryResults.txt');
        
        
        %% Load in summary results and parse
        
        % Column 1 is fUS computer timestamp
        % Column 2 is behavior computer timestamp
        % Column 3 is acquisition timestamp
        % Column 4 is phase
        % Column 5 is number of frames in memory
        % Column 6 is actual class
        % Column 7 is predicted class
        t_fUSacqs = textscan(fopen(summary_results_fname), '%s %s %s %d %d %d %d'); % load fUS acq time log - Most recent time log format
        fUS_comp_time = parse_fUS_time(t_fUSacqs{1});
        behavior_comp_time = parse_fUS_time(t_fUSacqs{2});
        fUS_acq_time = parse_fUS_time(t_fUSacqs{3});
        phase = double(cell2mat(t_fUSacqs(4)));
        memory_length = double(cell2mat(t_fUSacqs(5)));
        actual_class = double(cell2mat(t_fUSacqs(6)));
        predicted_class = double(cell2mat(t_fUSacqs(7)));
        
    else
        data = load_doppler_data('session_run_list', [session run]);
        phase = data.state_num;
        memory_length = data.memory_length;
        actual_class = data.actual_labels;
        predicted_class = data.predicted_labels;
        
        % In the data saved to MATLAB, zeros are represented by NaN.
        % Standardize by switching NaN back to zeros.
        phase(isnan(phase)) = 0;
        actual_class(isnan(actual_class)) = 0;
        predicted_class(isnan(predicted_class)) = 0;
    end
    
    
    %% Find number of training trials used (if any)
    if project_record.simulated_results(indx(session_counter))
        % For simulated runs, we can just load the experiment parameters
        % Load the experiment parameters
        experiment_parameters_fname = fullfile(base_path, subfolder_path, 'experiment_parameters.mat');
        exp_parameters = load(experiment_parameters_fname);
        number_training_trials(session_counter) = exp_parameters.nTrials_before_train;
        
        if project_record.pretraining(indx(session_counter))
            % If a previous session was used, then the training trials
            % were satisfied by the previous session.
            number_training_trials(session_counter) = 0;
        end
    else
        % This works if only successful trials are added to the training set.
        % If unsuccessful trials are also added, then this is not
        % perfect because it looks for certain types of trials. The
        % training set can be completed but the monkey still makes
        % behavioral mistakes. Therefore, this number is the upper
        % bound on # of training trials rather than the actual # of
        % training trials.
        
        % When was first prediction made
        first_prediction_ind = find(predicted_class, 1);
        
        % Find when the state changes
        state_change_ind = [diff(phase); 0]; % To keep same length as success_iti_ind
        
        % Where are the trials where predictions could have been made
        success_iti_ind = phase == 8;
        
        trial_ind = find(state_change_ind & success_iti_ind);
        number_training_trials(session_counter) = nnz(trial_ind < first_prediction_ind);
        
        % If we start part way through a trial, then may need to subtract one from
        % the training trials number. `bci_data_access` won't make prediction
        % unless it also sees a previous phase to memory (i.e. less than phase #4).
        
        first_new_trial_transition = find(state_change_ind < 0, 1);
        if phase(1) >= 4 && phase(first_new_trial_transition) == 8
            number_training_trials(session_counter) = number_training_trials(session_counter) - 1;
        end    
    end
    
    %% Analyze performance across time
    % Remove zeros
    zero_ind = predicted_class == 0;
    predicted_class_noZeros{session_counter} = predicted_class(~zero_ind);
    actual_class_noZeros{session_counter} = actual_class(~zero_ind);
    
    % Find number of correct choices
    correct_choices = predicted_class_noZeros{session_counter} == actual_class_noZeros{session_counter};
    nTrials = length(correct_choices);
    
    % Calculate performance across time
    performance_across_time = zeros(nTrials, 1);
    binomial_envelope = zeros(nTrials, 1);
    for i = 1:nTrials
        if isnan(rolling_window_size)
            performance_across_time(i) = nnz(correct_choices(1:i))/i;
            binomial_envelope(i) = (binoinv(1-pvalue_threshold, i, 1/nTargets) + 1)/i;
        else
            if i<=rolling_window_size
                start_trial = 1;
            else
                start_trial = i - rolling_window_size + 1;
            end
            end_trial = i;
            performance_across_time(i) = nnz(correct_choices(start_trial:end_trial))/min(i, rolling_window_size);
            binomial_envelope(i) = (binoinv(1-pvalue_threshold, min(i, rolling_window_size), 1/nTargets) + 1)/min(i, rolling_window_size);
        end
    end
    
    % If any above 1, then fix those.
    binomial_envelope(binomial_envelope > 1) = 1;
    
    % Find last time the decoder performance goes below desired pvalue
    % threshold.
    try
        last_nonsignificant_trial(session_counter) = find(performance_across_time-binomial_envelope<=0, 1, 'last');
    catch
        warning('No significant decoding');
        last_nonsignificant_trial(session_counter) = NaN;
    end
    
    % Create confusion matrix
    confusion_mat_end_trial = length(actual_class_noZeros{session_counter});
    if isnan(rolling_window_size)
        confusion_mat_start_trial = 1;
    else
        confusion_mat_start_trial = max(confusion_mat_end_trial-rolling_window_size+1, 1);
    end
    trial_range = confusion_mat_start_trial:confusion_mat_end_trial;
    if nTargets == 8
        % Multicoder creates an extra (center, center) group
        confusion_mat = build_confusion_matrix(nTargets+1, actual_class_noZeros{session_counter}(trial_range), predicted_class_noZeros{session_counter}(trial_range));
    else
        confusion_mat = build_confusion_matrix(nTargets, actual_class_noZeros{session_counter}(trial_range), predicted_class_noZeros{session_counter}(trial_range));
    end
    
    if update_confusion_matrix_for_each_trial
        figure(999);
        for trial_ind = 1:length(actual_class_noZeros{session_counter})
            switch nTargets
                case 2
                    confusion_mat_partial = build_confusion_matrix(nTargets, actual_class_noZeros{session_counter}(1:trial_ind), predicted_class_noZeros{session_counter}(1:trial_ind));
                    ph = create_confusion_matrix(confusion_mat_partial, ...
                        'label_strings', {'Left', 'Right'},...
                        'title', sprintf('Trial 1:%d - Confusion matrix', number_training_trials+trial_ind), ...
                        'FontColor', 'k');
                    drawnow;
                    pause(pause_length);
                    
                case 8
                    confusion_mat_partial = build_confusion_matrix(nTargets+1, actual_class_noZeros{session_counter}(1:trial_ind), predicted_class_noZeros{session_counter}(1:trial_ind));
                    confusion_chart = create_confusion_matrix(confusion_mat_partial([7 8 9 6 3 2 1 4], [7 8 9 6 3 2 1 4 5]), ...
                        'label_strings', reordered_label_strings, ...
                        'title', sprintf('Trial 1:%d - Confusion matrix', number_training_trials+trial_ind),...
                        'FontColor', 'k');
                    drawnow;
                    pause(pause_length);
            end
        end
    end
    
    
    confusion_mat_acrossSession(:, :, session_counter) = confusion_mat;
    performance_acrossSession{session_counter} = performance_across_time;
    
    % Find angular error for each session
    [empiric_mean_angular_error{session_counter}, null_mean_angular_error{session_counter}, pvalue_empiric{session_counter}] = calculate_angular_error_across_trials(predicted_class_noZeros{session_counter}, actual_class_noZeros{session_counter}, ...
        'num_replicates', num_replicates, ...
        'rolling_average', rolling_window_size);
    
    % Find last time the decoder performance goes below desired pvalue
    % threshold.
    try
        last_nonsignificant_trial_empiric(session_counter) = find(pvalue_empiric{session_counter}>=pvalue_threshold, 1, 'last');
    catch
        warning('No significant decoding');
        last_nonsignificant_trial_empiric(session_counter) = NaN;
    end
    
end

%% NaN-pad the entries. This makes it easy to compare across sessions.
% Only needs to be done if more than one session/run is being loaded
if num_sessions > 1
    len = max(cellfun('length',performance_acrossSession));
    b = cellfun(@(x)[x;NaN(len-size(x,1),1)],performance_acrossSession,'UniformOutput',false);
    performance_acrossSession_mat = [b{:}];
    
    len = max(cellfun(@(x) size(x, 2),null_mean_angular_error));
    b = cellfun(@(x)[x NaN(size(x, 1), len-size(x,2))],null_mean_angular_error,'UniformOutput',false);
    null_mean_angular_error_mat = cell2mat(b);
    
    len = max(cellfun(@(x) size(x, 2),empiric_mean_angular_error));
    b = cellfun(@(x)[x NaN(size(x, 1), len-size(x,2))],empiric_mean_angular_error,'UniformOutput',false);
    empiric_mean_angular_error_mat = cell2mat(b);
end

%% Calculate mean and ste
% Only performed if more than one session/run being loaded
% Since we are not ignoring NaN entries, this has the effect of reducing
% the analysis period to the session with the fewest trials.
if num_sessions > 1
    mean_performance = mean(performance_acrossSession_mat, 2);
    std_performance = std(performance_acrossSession_mat, 0, 2);
    ste_performance = std_performance/sqrt(size(performance_acrossSession, 2));
    
    
    % Remove NaNs
    nan_ind = isnan(mean_performance);
    mean_performance(nan_ind) = [];
    ste_performance(nan_ind) = [];
    
    % Calculate actual angular error
    mean_angular_error = mean(empiric_mean_angular_error_mat, 1);
    std_angular_error = std(empiric_mean_angular_error_mat, 0,1);
    ste_angular_error = std_performance/sqrt(size(empiric_mean_angular_error_mat, 1));
    
    
    % Remove NaNs
    nan_ind = isnan(mean_angular_error);
    mean_angular_error(nan_ind) = [];
    ste_angular_error(nan_ind) = [];
end

%% Calculate statistics for final performance
% Only performed if more than one session/run is being loaded
if num_sessions > 1
    [final_performance, final_angular_error] = deal(NaN(num_sessions, 1));
    for i = 1:num_sessions
        final_performance(i) = performance_acrossSession{i}(end);
        final_angular_error(i) = empiric_mean_angular_error{i}(end);
    end
    
    % Mean + SEM decoder accuracy
    performance_mean = mean(final_performance);
    performance_ste = std(final_performance)/sqrt(num_sessions);
    
    % Mean + SEM angular error
    angular_error_mean = mean(final_angular_error);
    angular_error_ste = std(final_angular_error)/sqrt(num_sessions);
    
    % How many trials did it take to reach significance
    trials_to_significance_mean = mean(last_nonsignificant_trial_empiric);
    trials_to_significance_ste = std(last_nonsignificant_trial_empiric)/sqrt(num_sessions);
    
    fprintf('Mean final accuracy is %0.2f +/- %0.2f%% \n', performance_mean*100, performance_ste*100);
    fprintf('Mean final angular error is %0.2f +/- %0.2f rad \n', angular_error_mean, angular_error_ste);
    fprintf('Mean trials to significance is %0.2f +/- %0.2f trials \n', trials_to_significance_mean, trials_to_significance_ste);
end



%% Plot performance across sessions
% Only performed if more than one session being loaded
if num_sessions > 1
    
    session_colors = distinguishable_colors(num_sessions);
    
    % Calculate last significant trial of the mean+ste plot
    nTrials = length(mean_performance);
    binomial_envelope = zeros(nTrials, 1);
    for i = 1:nTrials
        if isnan(rolling_window_size)
            binomial_envelope(i) = (binoinv(1-pvalue_threshold, i, 1/nTargets) + 1)/i;
        else
            if i<=rolling_window_size
                start_trial = 1;
            else
                start_trial = i - rolling_window_size + 1;
            end
            end_trial = i;
            binomial_envelope(i) = (binoinv(1-pvalue_threshold, min(i, rolling_window_size), 1/nTargets) + 1)/min(i, rolling_window_size);
        end
    end
    
    % If any above 1, then fix those.
    binomial_envelope(binomial_envelope > 1) = 1;
    
    % Find last time the decoder performance goes below desired pvalue
    % threshold.
    significance_threshold = ste_performance;
    last_nonsignificant_trial_meanSTE = find(mean_performance - significance_threshold - binomial_envelope<=0, 1, 'last');
    
    
    trial_numbers = 1:length(mean_performance);
    max_nTrials = 0;
    
    figure;
    tld = tiledlayout(1, 3);
    set(gcf,'color','k');
    nexttile();
    h1 = shadedErrorBar(trial_numbers, ...
        mean_performance,...
        ste_performance, ...
        'lineProps','b');
    ylabel('Percent correct');
    xlabel('Trial number (Aligned to first predicted trial)');
    title('Mean Performance Across Sessions ', 'Color', 'w');
    
    % Overlay the binomial envelope
    hold on;
    y = 1/nTargets * ones(nTrials, 1);
    chance_envelope = binomial_envelope - 1/nTargets;
    h2 = shadedErrorBar(1:length(y), y, chance_envelope);
    
    
    legend([h1.patch h2.patch], 'Decoder performance (Mean +/- STE)', sprintf('%0.2f%% binomial significance envelope', (1-pvalue_threshold)*100), 'Location', 'best');
    
    
    % Add text for chance level
    text(nTrials, 1/nTargets+0.015, ...
        'Chance performance', ...
        'Color', 'k', ...
        'HorizontalAlignment', 'right');
    
    if isempty(nTrials_for_plotting)
        xlim([1 nTrials]);
    else
        xlim(nTrials_for_plotting);
    end
    ylim([0 1]);
    set(gca, 'XColor', 'w');
    set(gca, 'YColor', 'w');
    title(tld, sprintf('Performance for Across Sessions for Monkey %s, Slot %s with %d Directions', monkey, slot, nTargets), ...
        'Color', 'w');
    
    
    nexttile();
    show_binomial_envelope = false;
    title('Performance in each session', 'Color', 'w');
    hold on;
    [h_performance, session_strings] = deal(cell(num_sessions, 1));
    h_significancemarker = cell(num_sessions, 1);
    for i = 1:num_sessions
        performance_across_time_plot = performance_acrossSession{i};
        nTrials = length(performance_across_time_plot);
        trial_numbers = number_training_trials(i)+1:number_training_trials(i)+length(performance_across_time_plot);
        max_nTrials = max([max_nTrials trial_numbers]);
        h_performance{i} = plot(trial_numbers, performance_across_time_plot, 'Color', session_colors(i, :));
        params_used = sprintf('Retrain: %d, PrevSess: %d', ...
            project_record.retraining(indx(session_counter)), ...
            project_record.pretraining(indx(session_counter)));
        session_strings{i} = sprintf('S%dR%d - %s -%s', session_run_list(i, 1), session_run_list(i, 2), params_used, datestr(project_record.Date(indx(session_counter)), 'YYYY-mmmm-dd'));
        
        if last_nonsignificant_trial(i) < nTrials+number_training_trials(i)
            h_significancemarker{i} = plot(number_training_trials(i)+last_nonsignificant_trial(i), performance_across_time_plot(last_nonsignificant_trial(i)), ...
                'MarkerEdgeColor', session_colors(i, :), 'Marker', 'd', ...
                'LineStyle', 'none');
        end
    end
    
    ylabel('Percent correct');
    xlabel('Trial number');
    chance_level = 1/nTargets * ones(max_nTrials, 1);
    if show_binomial_envelope
        
        binomial_envelope_maxTrials = zeros(max_nTrials, 1);
        for i = 1:max_nTrials
            binomial_envelope_maxTrials(i) = (binoinv(1-pvalue_threshold, i, 1/nTargets) + 1)/i;
        end
        
        % If any above 1, then fix those.
        binomial_envelope_maxTrials(binomial_envelope_maxTrials > 1) = 1;
        chance_envelope = binomial_envelope_maxTrials - 1/nTargets;
        h2 = shadedErrorBar(1:length(chance_level), chance_level, chance_envelope);
        
        legend([h_performance{:} h2.patch h_significancemarker{1}], session_strings{:}, sprintf('%0.2f%% binomial significance envelope', (1-pvalue_threshold)*100), 'Last nonsignificant trial per session', 'Location', 'best');
        
    else
        h2 = plot(1:length(chance_level), chance_level, 'k--');
        legend([h_performance{:} h_significancemarker{1}], session_strings{:}, 'Last nonsignificant trial per session', 'Location', 'best',...
            'Interpreter','none');
        
    end
    % Add text for chance level
    text(max_nTrials, 1/nTargets+0.015, ...
        'Chance performance', ...
        'Color', 'k', ...
        'HorizontalAlignment', 'right');
    
    
    if isempty(nTrials_for_plotting)
        xlim([1 max_nTrials]);
    else
        xlim(nTrials_for_plotting);
    end
    ylim([0 1]);
    set(gca, 'XColor', 'w');
    set(gca, 'YColor', 'w');
    title(tld, sprintf('Performance for Across Sessions for Monkey %s, Slot %s with %d Directions', monkey, slot, nTargets), ...
        'Color', 'w', ...
        'Interpreter', 'none');
    
    xtick_nums = xticks;
    xticks([1 xtick_nums]);
    
    % Angular error null distribution for each session on single plot
    nexttile()
    h_significancemarker_angular_error = cell(num_sessions, 1);
    
    
    hold on;
    h_actual_angular_error = cell(num_sessions, 1);
    for i = 1:num_sessions
        angular_mean_error = empiric_mean_angular_error{i};
        nTrials = length(angular_mean_error);
        
        h_actual_angular_error{i} = plot(1+number_training_trials(i):nTrials+number_training_trials(i), empiric_mean_angular_error{i},...
            'Color', session_colors(i, :));
        
        if last_nonsignificant_trial_empiric(i) < nTrials+number_training_trials(i)
            h_significancemarker_angular_error{i} = plot(number_training_trials(i)+last_nonsignificant_trial_empiric(i), angular_mean_error(last_nonsignificant_trial_empiric(i)), ...
                'MarkerEdgeColor', session_colors(i, :), 'Marker', 'd', ...
                'LineStyle', 'none');
        end
    end
    
    ylabel('Absolute mean angular error (rad)');
    xlabel('Trial number');
    set(gca, 'XColor', 'w');
    set(gca, 'YColor', 'w');
    yticks(0:pi/4:pi);
    yticklabels({'0', '\pi/4', '\pi/2', '3\pi/4', '\pi'});
    
    if isempty(nTrials_for_plotting)
        xlim([1 max_nTrials]);
    else
        xlim(nTrials_for_plotting);
    end
    
    ylim([0 pi]);
    
    title('Angular error in each session', 'Color', 'w');
    
    
    % Get null distribution
    % Need to retrieve this again rather than use the earlier null
    % distributions because the training trials may increase the necessary
    % number of trials for the null angular error to extend across entire
    % figure.
    [~, null_angular_error, ~] = calculate_angular_error_across_trials(randi(8, max_nTrials, 1), randi(8, max_nTrials, 1), ...
        'num_replicates', num_replicates);
    null_mean_angular_error = mean(null_angular_error);
    
    h_chance_empiric = plot(1:length(null_mean_angular_error), null_mean_angular_error, 'k--');
    legend([h_actual_angular_error{:} h_significancemarker_angular_error{1}], session_strings{:}, 'Last nonsignificant trial per session', 'Location', 'best', ...
        'Interpreter','none');
    % Add text for chance level
    text(max_nTrials, null_mean_angular_error(end)+0.015*pi, ...
        'Chance angular error', ...
        'Color', 'k', ...
        'HorizontalAlignment', 'right');
    
end


%% Plot each confusion matrix and across-trial performance

for i = 1:num_sessions
    
    actual_class = actual_class_noZeros{i};
    predicted_class = predicted_class_noZeros{i};
    
    subtitle = sprintf('S%dR%d', session_run_list(i, 1), session_run_list(i, 2));
    
    % Specify x limit for figures
    if isempty(nTrials_for_plotting)
        xlim_plot = [];
    else
        xlim_plot = nTrials_for_plotting;
    end
    
    % Specify ylimits for figures
    if isempty(accuracy_range_for_plotting)
        accuracy_ylim = [0 1];
    else
        accuracy_ylim = accuracy_range_for_plotting;
    end
    
    if isempty(angular_error_for_plotting)
        angular_error_ylim = [0 pi];
    else
        angular_error_ylim = angular_error_for_plotting;
    end
    
    % Specify colorbar limits for confusion matrix
    if isempty(accuracy_colorbar_limits)
        confusionmat_clim = [0 100];
    else
        confusionmat_clim = accuracy_colorbar_limits;
    end
    
    plot_performance(actual_class, predicted_class, ...
        'subtitle', subtitle, ...
        'pvalue_threshold', pvalue_threshold, ...
        'n_training_trials', number_training_trials(i), ...
        'num_replicates', num_replicates, ...
        'rolling_average', rolling_window_size, ...
        'xlim', xlim_plot, ...
        'accuracy_ylim', accuracy_ylim, ...
        'angular_error_ylim', angular_error_ylim, ...
        'confusionmat_clim', confusionmat_clim);
    
end


