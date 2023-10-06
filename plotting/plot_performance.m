function plot_performance(actual_class, predicted_class, varargin)
% plot_performance  Display consistent graphics for RT performance
%
% INPUTS:
%   actual_class:               (n x 1) double; The true labels
%   predicted_class:            (n x 1) double; The predicted labels
%   varargin: 
%       pvalue_threshold:       scalar double; Threshold for significance
%       n_training_trials:      scalar double; How many training trials?
%       subtitle:               string; Figure subtitle
%       num_replicates:         scalar double; How many replicates for
%                               angular error null distribution
%       rolling_average:        scalar; How many windows in the rolling
%                               average. NaN will use all avaiable trials
%                               in the average, i.e. no rolling average.
%       xlim:                   [x_min x_max]; x-axis limits
%       accuracy_ylim:          [y_min y_max]; y-axis limits for accuracy
%                               plot
%       angular_error_ylim:     [y_min y_max]; y-axis limits for angular
%                               error plot
%       confusionmat_clim:      [min max]; limits of colorbar for confusion
%                               matrix; between [0 100];
%      
% OUTPUTS:
%   None
%


%%
p = inputParser;
addOptional(p, 'pvalue_threshold', 0.05);
addOptional(p, 'n_training_trials', 0);
addOptional(p, 'subtitle', '');
addOptional(p, 'num_replicates', 1000);
addOptional(p, 'rolling_average', NaN); %Length of window
addOptional(p, 'xlim', []); % Trials to show
addOptional(p, 'accuracy_ylim', [0 1]); % Decimal
addOptional(p, 'angular_error_ylim', [0 pi]); %Radians
addOptional(p, 'confusionmat_clim', [0 100]); % Percentage
parse(p, varargin{:});
inputs = p.Results;

%% Find how many classes
nClasses = length(unique(actual_class));

%% Analyze performance across time
pvalue_threshold = inputs.pvalue_threshold;


% Remove zeros
zero_ind = predicted_class == 0;
predicted_class_clean = predicted_class(~zero_ind);
actual_class_clean = actual_class(~zero_ind);

correct_choices = predicted_class_clean == actual_class_clean;
nTrials = length(correct_choices);

performance_across_time = zeros(nTrials, 1);
binomial_envelope = zeros(nTrials, 1);
for i = 1:nTrials
    if isnan(inputs.rolling_average)
        performance_across_time(i) = nnz(correct_choices(1:i))/i;
        binomial_envelope(i) = (binoinv(1-pvalue_threshold, i, 1/nClasses) + 1)/i;
    else
        if i<=inputs.rolling_average
            start_trial = 1;
        else
            start_trial = i - inputs.rolling_average + 1;
        end
        end_trial = i;
        performance_across_time(i) = nnz(correct_choices(start_trial:end_trial))/min(i, inputs.rolling_average);
        binomial_envelope(i) = (binoinv(1-pvalue_threshold, min(i, inputs.rolling_average), 1/nClasses) + 1)/min(i, inputs.rolling_average);
    end
    

end

% If any above 1, then fix those.
binomial_envelope(binomial_envelope > 1) = 1;

% Find last time the decoder performance goes below desired pvalue
% threshold.
last_nonsignificant_trial = find(performance_across_time-binomial_envelope<=0, 1, 'last');

% Create confusion matrix

confusion_mat_end_trial = length(actual_class_clean);
if isnan(inputs.rolling_average)
   confusion_mat_start_trial = 1;
else
   confusion_mat_start_trial = max(confusion_mat_end_trial-inputs.rolling_average+1, 1);
end
trial_range = confusion_mat_start_trial:confusion_mat_end_trial;

if nClasses == 8
    % Multicoder creates an extra (center, center) group
    confusion_mat = zeros(nClasses+1, nClasses+1);
    for ac = 1:nClasses+1
        for pr = 1:nClasses+1
            confusion_mat(ac, pr) = nnz(actual_class_clean(trial_range) == ac & predicted_class_clean(trial_range) == pr);
        end
    end
else
    confusion_mat = zeros(nClasses, nClasses);
    for ac = 1:nClasses
        for pr = 1:nClasses
            confusion_mat(ac, pr) = nnz(actual_class_clean(trial_range) == ac & predicted_class_clean(trial_range) == pr);
        end
    end
end

%% Plotting

number_training_trials = inputs.n_training_trials;
figure('Units', 'normalized', 'Position', [0.15 0.25 0.7 0.5]);
tld = tiledlayout(1,3);
nexttile()
h1 = plot(number_training_trials+1:number_training_trials+length(performance_across_time), performance_across_time);
ylabel('Percent correct');
xlabel('Trial number');
set(gca, 'XColor','w', 'YColor','w')


% Overlay the binomial envelope
hold on;
y = 1/nClasses * ones(nTrials, 1);
envelope_y = binomial_envelope - 1/nClasses;
h2 = shadedErrorBar(number_training_trials+1:number_training_trials+length(y), y, envelope_y);


if last_nonsignificant_trial < nTrials
    % Mark the last nonsignificant trial
    h3 = plot([number_training_trials+last_nonsignificant_trial number_training_trials+last_nonsignificant_trial], [0 1], 'r');
    
    legend([h1 h2.patch h3], 'Decoder performance', sprintf('%0.2f%% binomial significance envelope', (1-pvalue_threshold)*100), 'Last nonsignificant trial', 'Location', 'best');
    
    % Add text
    text(number_training_trials+last_nonsignificant_trial + nTrials*0.01, 0.1, ...
        sprintf('Trial %d', number_training_trials+last_nonsignificant_trial), ...
        'Color', 'r');
else
    legend([h1 h2.patch], 'Decoder performance', sprintf('%0.2f%% binomial significance envelope', (1-pvalue_threshold)*100), 'Location', 'best');
end


if isempty(inputs.xlim)
    xlim([1+number_training_trials nTrials+number_training_trials]);
else
    xlim(inputs.xlim);
end

ylim(inputs.accuracy_ylim);

%% Plot confusion matrix
nexttile();

switch nClasses
    case 2
        confusion_chart = create_confusion_matrix(confusion_mat,...
            'label_strings', {'Left', 'Right'},...
            'title', 'Combined prediction', ...
            'FontColor', 'w');
        set(gcf,'color','k');
        axis image;
    case 8
            combined_label_strings = {'Ipsilateral Down', 'Center Down', 'Contralateral Down', ...
        'Ipsilateral Center',  'Middle', 'Contralateral Center', ...
        'Ipsilateral Up', 'Center Up', 'Contralateral Up'};
        reordered_label_strings = combined_label_strings([7 8 9 6 3 2 1 4 5]);
        
        reordered_confusion_mat = confusion_mat([7 8 9 6 3 2 1 4 5], [7 8 9 6 3 2 1 4 5]);

        confusion_chart = create_confusion_matrix(reordered_confusion_mat(1:8, :), ...
            'label_strings', reordered_label_strings, ...
            'title', 'Combined prediction', ...
            'FontColor', 'w', ...
            'colorbar_limits', inputs.confusionmat_clim);
        set(gcf,'color','k');
        axis image;
        
end

%% Calculate mean angular error and null distribution of angular errors
[empiric_mean_angular_error, null_mean_angular_error, pvalue_empiric] = calculate_angular_error_across_trials(predicted_class_clean, actual_class_clean, ...
    'num_replicates', inputs.num_replicates, ...
    'rolling_average', inputs.rolling_average);

mean_null_distribution_across_trials = mean(null_mean_angular_error, 1);
max_nTrials = size(mean_null_distribution_across_trials, 2);
if ~isnan(inputs.rolling_average)
    mean_null_distribution_across_trials(inputs.rolling_average+1:end) = mean_null_distribution_across_trials(min(inputs.rolling_average, max_nTrials));
end
quantiles_null_distribution_across_trials = quantile(null_mean_angular_error, [pvalue_threshold 1-pvalue_threshold]);
if ~isnan(inputs.rolling_average)
    quantiles_null_distribution_across_trials(:, inputs.rolling_average+1:end) = repmat(quantiles_null_distribution_across_trials(:, min(inputs.rolling_average, max_nTrials)), ...
        1, ...
        size(quantiles_null_distribution_across_trials(:, inputs.rolling_average+1:end), 2));
end

error_bars = abs(quantiles_null_distribution_across_trials([2 1],:) - mean_null_distribution_across_trials);

nexttile();
hold on;
h_null_distribution = shadedErrorBar(number_training_trials+1:number_training_trials+nTrials, mean_null_distribution_across_trials , error_bars);
h_actual_angular_error = plot(number_training_trials+1:number_training_trials+nTrials, empiric_mean_angular_error);
ylabel('Mean absolute angular error (deg)');
xlabel('Trial number');
yticks(0:pi/4:pi);
yticklabels({'0', '45', '90', '135', '180'});

if isempty(inputs.xlim)
    xlim([1+number_training_trials nTrials+number_training_trials]);
else
    xlim(inputs.xlim);
end

ylim(inputs.angular_error_ylim);


legend([h_actual_angular_error h_null_distribution.patch], 'Decoder performance', sprintf('%0.2f%% significance envelope', (1-pvalue_threshold)*100), 'Location', 'best');

set(gca, 'XColor','w', 'YColor','w')


%% Add overall title
title(tld, ['Decoder performance - ' inputs.subtitle], 'Color', 'w');

