function numTime = parse_fUS_time(chartime)
% parse_fUS_time  Converts the character arrays in the behavior data to
% numbers 

% INPUTS:
%   chartime:       (n x 1) cell array where each cell is the exact same 
%                   length char array
%      
% OUTPUTS:
%   numTime:        (n x 1) double; time in seconds.
%
% See also parse_behavior_time

%% Parse fUS time
% Since each cell entry is the exact same length, we can convert to a more
% convenient format 
chartime = cell2mat(chartime);

% extract hours/min/sec/ms from chartimes
HH = str2num(chartime(:,1:2));   % Hours array
MM = str2num(chartime(:,4:5));   % Minutes array
SS = str2num(chartime(:,7:8));   % seconds array
MS = str2num(chartime(:,10:12)); % milliseconds array

% concat visual stim trial times [HH MM SS.MS]
numTimeArray = cat(2, HH, MM, SS+1e-3*MS);

% convert time to [s]
numTime =  numTimeArray(:,1)*3600 + numTimeArray(:,2)*60+ numTimeArray(:,3);
end