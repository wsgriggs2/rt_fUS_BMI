function path = get_data_path(varargin)
% get_data_path Retrieve path to the user's data
%
% INPUTS:
%   varargin:
%       path_type:      String; Available options are: 
%                       {'root', 'doppler', 'sulcus', 'output', 
%                       'simulated','aligned'}
% 
% Outputs:
%   path:               string; Absolute path to where fUS data is stored

%% Input parser
p = inputParser;
p.addParameter('path_type','default', @ischar)
p.parse(varargin{:});
inputs = p.Results;

%% Check if root has already been set
if isfile('data_path.mat')
    % If path has already been specified, then load that.
    load('data_path.mat', 'base_path');
else
    % If path has not been set, then ask user to specify a location
    base_path = specify_data_path;
end

%% create full path 
switch inputs.path_type
    case {'default', 'root'}
        path_suffix = '';
        
    case 'doppler'
        path_suffix = 'doppler data';
        
    case 'sulcus'
        path_suffix = 'sulcus';
        
    case 'output'
        path_suffix = 'output';
        
    case 'simulated'
        path_suffix = 'simulated';
        
    case 'aligned'
        path_suffix = 'aligned doppler data';
        
    otherwise
        error('Did not recognize path type');
end

% Return complete path
path = fullfile(base_path, path_suffix);


end