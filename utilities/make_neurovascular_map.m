function neurovascular_map = make_neurovascular_map(DopplerImages, varargin)
% make_neurovascular_map Generate a neurovascular image from raw Power 
% Doppler images
%
% INPUTS:
%   DopplerImage:           3D or 4D matrix; Power Doppler images
%                           3D (zPix x xPix x timepoints)
%                           4D (zPix x xPix x timepoints x trials)
%
%   varargin:
%       root:               double; What power of root do you want to apply
%                           to image?
%
% Outputs:
%   neurovascular_map:      2D matrix; Average neurovascular activity in 
%                           image


%% Variable inputs
p = inputParser;
p.addOptional('root',3)
p.parse(varargin{:});
inputs = p.Results;

%% Make the neurovascular map

if ndims(DopplerImages) == 4
    % If 4D array, then take average over timepoints, then across trials.
    neurovascular_map = nthroot(squeeze(mean(mean(DopplerImages,3),4)),inputs.root);

elseif ndims(DopplerImages) == 3
    % If 3D array, then take average across all timepoints.
    neurovascular_map = nthroot(squeeze(mean(DopplerImages, 3)), inputs.root);
elseif ismatrix(DopplerImages)
    % If just a single image, then apply the nth root to scale the
    % intensity values.
    neurovascular_map = nthroot(DopplerImages, inputs.root);
end

% Remove the bottom 1% of voxels to help with dynamic range.
lowerCutoff = quantile(neurovascular_map(:),0.01);
neurovascular_map(neurovascular_map<lowerCutoff) = lowerCutoff;
    
end