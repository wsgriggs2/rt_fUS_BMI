function alignment_tform = align_neurovascular_maps(neurovascular_map1, neurovascular_map2, varargin)
% align_neurovascular_maps  Find the rigid-body transform from map2 to map1
%
% INPUTS:
%   neurovascular_map1 -    2D matrix; Neurovascular map #1
%   neurovascular_map2 -    2D matrix; Neurovascular map #2
%
% OUTPUTS:
%   alignment_tform -       3x3 matrix; 2D affine transorm matrix -
%                           See https://en.wikipedia.org/wiki/Affine_transformation
%                           for description of 2D affine transform matrix
%

%% Variable inputs
p = inputParser;
addParameter(p, 'neurovascular_map1_session_run', []);
addParameter(p, 'neurovascular_map2_session_run', []);
parse(p, varargin{:});
inputs = p.Results;

%% Perform alignment
% Z-score each image to keep to the same range of values. Helps avoid the
% automated alignment algorithm from diverging.
neurovascular_map1_zscore = zscore(neurovascular_map1, 0, 'all');
neurovascular_map2_zscore = zscore(neurovascular_map2, 0, 'all');

% Align neurovascular_map2 to neurovascular_map1
app = image_alignment_app(neurovascular_map1_zscore, neurovascular_map2_zscore, ...
    'image1_session_run', inputs.neurovascular_map1_session_run, ...
    'image2_session_run', inputs.neurovascular_map2_session_run);
pauseCount=0;

% Mandatory wait period because other the code
% keeps running forward without waiting for user input
% from the manual image alignment GUI.
while isempty(app.alignment_tform)
    pause(1);
    pauseCount=pauseCount+1;
    if pauseCount == 2000
        error('Alignment_tform was not created within 2000 seconds. Stopping code');
        break;
    end
end
alignment_tform = app.alignment_tform;
close(app.UIFigure);
clear app;
