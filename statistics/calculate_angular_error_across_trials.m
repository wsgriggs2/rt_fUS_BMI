function [empiric_mean_angular_error, null_mean_angular_error, pvalue, empiric_angular_error, null_angular_error] = calculate_angular_error_across_trials(predicted_class, actual_class, varargin)
% calculate_angular_error_across_trials  Generates mean angular error 
% across trials and also generates null distribution for angular error
%
% INPUTS:
%   predicted_class:              (nx1); The predicted class labels
%   actual_class:                 (nx1); The actual class labels
%   varargin
%       num_replicates:           (scalar); How many replicates for the null
%                                 distribution
%       rolling_average:          scolar double; How many trials to perform
%                                 rolling average across. If NaN, then use
%                                 all avaiable trials, i.e. no rolling 
%                                 average.
%
% OUTPUTS:
%   empiric_mean_angular_error:   Actual mean angular error seen in
%                                 experiment
%   null_mean_angular_error:      Null distribution for angular error
%                                 assuming a uniform distribution.
%   pvalue:                       Empiric p-value. What percentage of the
%                                 null distribution exceeds the empiric 
%                                 mean angular error.
%   empiric_angular_error:        Angular error for a given trial
%   null_angular_error:           Angular error for each trial and
%                                 replicate


%% Handling varargin
p = inputParser;
p.addOptional('num_replicates',1000);
p.addOptional('rolling_average', NaN);
p.parse(varargin{:});
inputs = p.Results;


%% Clean up data depending on nClasses
% Find how many classes there are.
nClasses = length(unique(actual_class));

switch nClasses
    case 8
        % Find all the middle predictions
        middle_ind = predicted_class==5;
        

        % Assign middle predictions to the worst possible angular error
        % (pi).
        direction_predicted = predicted_class;
        direction_actual = actual_class;

        % Map back to 1:8 instead of 1:9
        direction_predicted(direction_predicted>5) = direction_predicted(direction_predicted > 5) - 1;
        direction_actual(direction_actual > 5) = direction_actual(direction_actual > 5) - 1;

        % Assign the middle ind to actual+pi for prediction. This is
        % worst case angular error.
        direction_predicted(middle_ind) = mod(direction_actual(middle_ind) + 4, 8);
        direction_predicted(direction_predicted == 0) = 8;

        
        % Map from predicted number to direction
        angle_lookup_table = 0:pi/4:7*pi/4;
        angle_lookup_table = angle_lookup_table([6 7 8 5 1 4 3 2]);
        
    case 2
        direction_predicted = predicted_class;
        direction_actual = actual_class;
        
        % Map from predicted number to direction
        angle_lookup_table = [0 pi];
        
end

%% Load null distribution
[null_angular_error, null_mean_angular_error] = load_null_angular_error_distribution(nClasses, inputs.num_replicates);

nTrials = length(direction_actual);

% Should already be loaded, we only need the first nTrials

null_mean_angular_error = null_mean_angular_error(:, 1:nTrials);
null_angular_error = null_angular_error(:, 1:nTrials);


% Pre-allocate space
[empiric_mean_angular_error, pvalue] = deal(NaN(1, nTrials));

% Calculate angular error for each trial
empiric_angular_error = abs(angdiff(angle_lookup_table(direction_predicted), angle_lookup_table(direction_actual)));
for trial = 1:nTrials
    % Calculate angular error and mean angular_error for the actual data
    % Probably a more efficient way to do this that avoids a for loop
    if isnan(inputs.rolling_average)
        empiric_mean_angular_error(trial) = mean(empiric_angular_error(1:trial));

        % Calculate percentage of null distribution exceeding true value
        pvalue(trial) = 1 - nnz(null_mean_angular_error(:, trial) > empiric_mean_angular_error(trial))/inputs.num_replicates;
    else
        if trial<=inputs.rolling_average
            start_trial = 1;
        else
            start_trial = trial - inputs.rolling_average + 1;
        end
        end_trial = trial;
        empiric_mean_angular_error(trial) = mean(empiric_angular_error(start_trial:end_trial));

        % Calculate percentage of null distribution exceeding true value
        pvalue(trial) = 1 - nnz(null_mean_angular_error(:, min(trial, inputs.rolling_average)) > empiric_mean_angular_error(trial))/inputs.num_replicates;
    end
end



