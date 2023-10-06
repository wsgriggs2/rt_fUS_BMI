function project_record = get_project_record(varargin)
% get_project_record Load and display project record
% 
% INPUTS:
%   varargin:
%       verbose:                bool; print project record to screen?
%       include_simulated:      bool; include simulated results in the
%                               project_record returned?
%
% OUTPUTS:
%   project_record:             (n x m) table; Each row represents
%                               information associated with one
%                               session/run, i.e. one experimental run.

p = inputParser;
p.addOptional('verbose', false);
p.addOptional('include_simulated', false);
p.parse(varargin{:});
inputs = p.Results;

filepath = get_data_path('path_type', 'root');
filename = 'project_record.json';


full_filename = fullfile(filepath, filename);

if ~isfile(full_filename)
    initial_loc = fullfile(filepath, '*.json');
    title_string = ['Could not automatically find project record.'...
        'Please manually select Project Record JSON file'];
    disp(title_string);
    [filename, filepath] = uigetfile(initial_loc, ...
        title_string);
    full_filename = fullfile(filepath, filename);
end


project_record = load_json_as_table(full_filename);

if ~inputs.include_simulated
    simulated_indx = project_record.simulated_results;
    project_record = project_record(~simulated_indx, :);
end

if inputs.verbose
    disp(project_record)
end