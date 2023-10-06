function [null_angular_error, null_mean_angular_error] = generate_null_distribution_angular_error(num_classes, num_replicates)
% generate_null_distribution_angular_error  Generate the null distribution
% for angular error by randomly sampling from a uniform distribution of the
% possible directions.
%
% INPUTS:
%   num_classes:                double; 2 or 8
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
% See also load_null_angular_error_distribution and
% load_null_mean_angualr_error

%% Set defaults
if ~exist('num_classes', 'var')
    num_classes = [2 8];
end

if ~exist('num_replicates', 'var')
    num_replicates = 1000;
end

%% Generate null distribution
for i = 1:length(num_classes)
    switch num_classes(i)
        case 8
            % Map from predicted number to direction
            angle_lookup_table = 0:pi/4:7*pi/4;
            angle_lookup_table = angle_lookup_table([6 7 8 5 1 4 3 2]);
        case 2
            % Map from predicted number to direction
            angle_lookup_table = [0 pi];
        otherwise
            error('This number of classes is not supported.');
    end
    
    
    num_trials = 2000; % Upper bound of number of trials ever expected
    
    % Pre-allocate space
    null_mean_angular_error = NaN(num_replicates, num_trials);
    
    
    % Calculate angular error for each trial
    % Randomly sample form uniform distribution
    random_direction = randi(num_classes(i), num_replicates, num_trials);
    null_angular_error = abs(angdiff(angle_lookup_table(random_direction), zeros(num_replicates, num_trials)));
    
    ppm = ParforProgressbar(num_trials*num_replicates);
    parfor trial = 1:num_trials
        % Calculate iterative mean as we add trials.
        for rep = 1:num_replicates
            null_mean_angular_error(rep, trial) = mean(null_angular_error(rep, 1:trial));
            ppm.increment();
        end
    end
    delete(ppm);
    
    %% Save the null distribution
    save(sprintf('angular_error_null_distribution_%dtgts_%dreps.mat', num_classes(i), num_replicates), 'null_angular_error', 'null_mean_angular_error');
    
end