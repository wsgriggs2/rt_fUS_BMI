function null_mean_angular_error = load_null_mean_angular_error(n_targets, num_replicates, n_trials)
% LOAD_NULL_MEAN_ANGULAR_ERROR  Helper function to load null mean angular 
% error. If requested variable not available, then ask user if they want to
% create the requested null distribution.
%
% INPUTS:
%   n_targets:                  scalar double;
%   num_replicates:             scalar double;
%   n_trials:                   scalar double;
%      
% OUTPUTS:
%   null_mean_angular_error:    (num_replicates x n_trials) double;
%                               Trial-wise iterative mean of the randomly 
%                               angular error.
%
% See also LOAD_NULL_ANGULAR_ERROR_DISTRIBUTION
% GENERATE_NULL_DISTRIBUTION_ANGULAR_ERROR

try
    null_distribution = load(sprintf('angular_error_null_distribution_%dtgts_%dreps.mat', n_targets, num_replicates), 'null_angular_error', 'null_mean_angular_error');
    null_mean_angular_error = null_distribution.null_mean_angular_error(:, 1:n_trials);
catch
    choice = questdlg('This combo of nClasses and num_replicates did not exist. Do you want to create now?',...
        'Create new null distribution?','yes','no','yes');
    if strcmp(choice,'yes')
        [~, null_mean_angular_error] = generate_null_distribution_angular_error(n_targets, num_replicates);
    else
        error('The distribution does not exist currently. Aborting.\n');
    end
    null_mean_angular_error = null_mean_angular_error(:, 1:n_trials);
end
              