function numTime = parse_behavior_time(behavior_struct, columnTitle)
% parse_behavior_time  Converts the character arrays in the behavior data
% to numbers
%
% INPUTS:
%   behavior_struct:    (1 x n) struct; Behavioral record
%   columnTitle:        string; Name of the behavioral state to convert.
%
% OUTPUTS:
%   numTime:            (n x 1) double; Time in seconds


%% Get the specified column
temp_cell = {behavior_struct.(columnTitle)}';
empty_ind = cellfun(@isempty, temp_cell);
negative_ind = false(size(empty_ind));
negative_ind(~empty_ind) = contains(temp_cell(~empty_ind), '-');

number_chars = size(temp_cell{1}, 2);
if any(negative_ind)
    chartime(~empty_ind & ~negative_ind, 2:number_chars) = cell2mat(temp_cell(~empty_ind & ~negative_ind));
    chartime(~empty_ind & negative_ind, :) = cell2mat(temp_cell(~empty_ind & negative_ind));
else
    chartime = cell2mat(temp_cell(~empty_ind));
end



% For flexibility, we are parsing this a single time string at a time.
numTimeArray = zeros(size(chartime, 1), 3);

for i = 1:size(chartime, 1)
    temp_chartime = string(chartime(i, :));
    colon_ind = strfind(temp_chartime, ':');
    comma_ind = strfind(temp_chartime, ',');
    negative_ind(i) = contains(temp_chartime, '-');
    numTimeArray(i, 1) = str2num(chartime(i,1:colon_ind(1)-1));   % Hours
    numTimeArray(i, 2) = str2num(chartime(i,colon_ind(1)+1:colon_ind(2)-1));   % Minutes
    numTimeArray(i, 3) = str2num(chartime(i,colon_ind(2)+1:comma_ind-1));   % seconds
    numTimeArray(i, 3) = numTimeArray(i, 3) + 1e-3*str2num(chartime(i,comma_ind+1:comma_ind+3)); % milliseconds
end


% Keep same size as input
numTime = NaN(max(size(behavior_struct)), 1);

% convert time to [s]
numTime(~empty_ind) =   numTimeArray(:,1)*3600 + numTimeArray(:,2)*60+ numTimeArray(:,3);
numTime(negative_ind) = numTime(negative_ind)*-1;
end