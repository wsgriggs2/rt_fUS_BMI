function [null_angular_error, null_mean_angular_error] = load_null_angular_error_distribution(num_classes, num_replicates)
% LOAD_NULL_ANGULAR_ERROR_DISTRIBUTION  Load the null distributions from 
% file. If the requested null distribution does not exist yet, then ask if
% the user wants to create it.
%
% INPUTS:
%   num_classes:                scalar double; 2 or 8
%   num_replicates:             scalar double; (n); How many replicates?
%      
% OUTPUTS:
%   null_angular_error:         (n x m) double; Randomly sampled angular
%                               error drawn from uniform distribution for
%                               each trial. n = num_replicates. m =
%                               num_trials
%   null_mean_angular_error:    (n x m) double; Trial-wise iterative mean
%                               of the randomly angular error. 
%                               n = num_replicates. m = num_trials
%
% See also LOAD_NULL_MEAN_ANGULAR_ERROR,
% GENERATE_NULL_DISTRIBUTION_ANGULAR_ERROR


try
    load(sprintf('angular_error_null_distribution_%dtgts_%dreps.mat', num_classes, num_replicates));
catch
    fprintf('The angular error null distribution for %d targets does not exist yet.\n', num_classes);
    fprintf('Generating null distribution now with %d replicates.\n', num_replicates);
    [null_angular_error, null_mean_angular_error] = generate_null_distribution_angular_error(num_classes, num_replicates);
    
    
end
end