function add_entry_to_project_record(new_entry)
% add_entry_to_project_record  Add entry to `project_record.json`. If the
% JSON file does not exist, then create it.
%
% INPUTS:
%   new_entry:       (1 x 1) table; New entry for project_record.json. Must
%                    match existing table fields.
%      
% OUTPUTS:
%   None


%% Find automatically
filepath = get_data_path('path_type', 'root');
filename = 'project_record.json';

full_filename = fullfile(filepath, filename);

if ~isfile(full_filename)
    create_new_file = true;
else
    create_new_file = false;
end
%%
if create_new_file
    project_record = new_entry;
else
    project_record = load_json_as_table(full_filename);

    project_record = [project_record; new_entry];

    % Sort by session number and run
    project_record = sortrows(project_record, [1, 2]);
end

% Save table to file
save_table_as_json(project_record, fullfile(full_filename));
fprintf('DONE! %s has been updated\n\n\n', filename)