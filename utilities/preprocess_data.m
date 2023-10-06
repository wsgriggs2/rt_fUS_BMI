function dop_out = preprocess_data(dop_in, varargin)
% preprocess_data  Perform basic pre-processing steps on the doppler
% data of the form: [yPixels x xPixels x timeWindows]
%
% INPUTS:
%   dop_in:                 (yPixels x xPixels x nImages); Doppler data
%   varargin: 
%       zscore:             bool; Apply z-scoring?
%       spatial_filter:     cell array; Parameters that will be eventually
%                           passed to FSPECIAL.
%       verbose:            bool; Display output?
%      
% OUTPUTS:
%   dop_out:                Same size as dop_ins

%% handling varargin
p = inputParser;
p.addOptional('zscore',false,@islogical)
p.addOptional('spatial_filter',{'disk', 2, 0})
p.addOptional('verbose',false,@islogical)
p.parse(varargin{:});
sets = p.Results;
if sets.verbose
    disp(sets)
end

%% main function is here (function calls)
if sets.zscore, dop_in = dopZ(dop_in); end
if any(~cellfun(@isempty, sets.spatial_filter)), dop_in = spatial_filter(dop_in, sets.spatial_filter); end

dop_out = dop_in;

%% spatial filter
    function data_out = spatial_filter(data_in, filter_options)
        % Apply specified spatial filter
        % First entry in `filter_options` cell defines the spatial filter
        % Subsequent entries define the parameters for that specified
        % filter.
        %
        % For 'disk' - specify filter size
        % For 'gaussian - specify filter size and filter sigma
        [~, ~, n_windows] = size(data_in);
        data_out = NaN(size(data_in));
        switch filter_options{1}
            case 'disk'
                h = fspecial(filter_options{1}, filter_options{2}); % disk size defined here
            case 'gaussian'
                h = fspecial(filter_options{1}, filter_options{2}, filter_options{3}); % disk size defined here
            otherwise
                error('This filter type has not been implemented');
        end
        
        
        % can't filter2 n-d arrays -> ugly nested for-loop (sorry)
        % note: NOT faster with parfor overhead (I tried)
        
        for window = 1:n_windows
            data_out(:, :, window) = ...
                filter2(h, squeeze(data_in(:, :, window)));
            
        end
    end


%% z score each voxel
    function data_out = dopZ(data_in)
        [n_depth, n_width, n_windows] = size(data_in);
        
        % flatten the data into a voxel-column for zscoring
        data_out = permute(data_in, [2, 1, 3]);% [n_width, n_depth, n_trials, n_windows]);
        data_out = reshape(data_out, [n_depth*n_width, n_windows]);
        % perform the zscore
        data_out = zscore(data_out')';
        % reshape back to original data dims
        data_out = reshape(data_out, [n_width, n_depth, n_windows]);
        data_out = permute(data_out, [2, 1, 3]);
    end


end