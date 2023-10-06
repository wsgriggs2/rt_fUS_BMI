function [behavior, timestamps] = clean_behavior_and_timestamps(behavior, timestamps)
% CLEAN_BEHAVIOR_AND_TIMESTAMPS  Perform preprocessing on behavioral
% record. Ensure that behavior is sorted by increasing trial number. Adjust
% behavior and timestamps so that they start from 0 and are directly
% comparable
%
% INPUTS:
%   behavior:           (1 x n) struct; The behavior record where each
%                       entry is one trial. The timestamps in this record
%                       are generally by default in time of day rather than
%                       time since start of session.
%   timestamps:         (n x 1) cell array; Time of each fUS frame
%
% OUTPUTS:
%   behavior:           The adjusted behavior record
%   timestamps:         The adjusted timestamps
%

%% Ensure behavior is sorted
[~, sort_ind] = sort([behavior.trialID]);
if any(diff(sort_ind) < 0)
    behavior = behavior(sort_ind);
end

%% Adjust timing in behavior, timestamps

 state_names = {'trialstart', 'initialfixation', 'fixationhold', 'cue', ...
        'memory', 'target_acquire', 'target_hold', 'reward', 'iti'}; 

minTime = min([parse_behavior_time(behavior(1), 'trialstart'), timestamps(:, 1)'], [], 'all');
timestamps(:, 1) = timestamps(:, 1) - minTime;
behavior = adjust_behavior_time(behavior, state_names, -minTime);