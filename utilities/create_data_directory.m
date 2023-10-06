function create_data_directory()
% create_data_directory Set up desired data folder organization
% This function will create the expected file directory organization

% Types of folders to create
folder_names = {'doppler', 'sulcus', 'output', 'simulated', 'aligned'};

%% Iterate over folder types and create folder if necessary
for folder_ind = 1:length(folder_names)
    pathname = get_data_path('path_type', folder_names{folder_ind});
    if ~isfolder(pathname)
        mkdir(pathname);
    end
end

%% If necessary, copy downloaded Doppler data to appropriate location
% Ask user if they want to copy data they already downloaded
% ask the user if they want to save

% Using copy instead of move to avoid any potentially destructive
% operations.
text_str = '--------------\nWould you like to copy previously downloaded Doppler data?\nyes/no: ';
copy_data = lower(input(text_str, 's'));
acceptable_inputs = {'yes', 'no'};
while ~any(strcmp(copy_data, acceptable_inputs))
    warning('Input must be "yes" or "no". Please try again.')
    copy_data = lower(input(text_str, 's'));
end

if strcmp(copy_data, 'yes')
    %Specify location of data
    title_string = 'Specify location of data to be copied';
    fprintf(['\n--------------\n' title_string '\n']);
    downloaded_data_loc = uigetdir(matlabroot, title_string);

    % Copy doppler data to new location
    disp('Copying Doppler data:');
    rt_fUS_data_file_info = dir(fullfile(downloaded_data_loc, 'rt_fUS_data_S*_R*.mat'));
    rt_fUS_data_file_info = natsortfiles(rt_fUS_data_file_info);
    new_folder_loc = get_data_path('path_type', 'doppler');
    for file_ind = 1:length(rt_fUS_data_file_info)
        filename =  rt_fUS_data_file_info(file_ind).name;
        fprintf('\n%s', filename);
        copyfile(fullfile(rt_fUS_data_file_info(file_ind).folder,  filename),...
            fullfile(new_folder_loc, filename));
    end
    disp('\nFinished copying Doppler data!');

    % Specify location of session record
    title_string = 'Select the downloaded Project Record JSON file';
    fprintf(['\n--------------\n' title_string '\n']);
    [project_record_filename, project_record_loc] = uigetfile(...
        fullfile(downloaded_data_loc, 'project_record.json'), ...
        title_string);

    % Specify location of data description pdf
    title_string = 'Select the downloaded file description PDF';
    fprintf(['\n--------------\n' title_string '\n']);
    [data_descriptor_filename, project_record_loc] = uigetfile(...
        fullfile(downloaded_data_loc, 'DescriptionOfFiles.pdf'), ...
        title_string);

    % Copy session record and description of data variables to new location
    disp('Copying Project Record JSON file');
    copyfile(fullfile(project_record_loc, project_record_filename), ...
        fullfile(get_data_path('path_type', 'root'), project_record_filename));
    disp('Finished copying Project Record JSON file!');

    disp('Copying data descriptor PDF file');
    copyfile(fullfile(project_record_loc, data_descriptor_filename), ...
        fullfile(get_data_path('path_type', 'root'), data_descriptor_filename));
    disp('Finished copying data descriptor PDF file!');

    % Ask if user wants to delete the original copies of the files now that
    % everything has been copied successfully
    text_str = '--------------\nWould you like to delete the original copy of the data you downloaded?\nyes/no: ';
    delete_orig_data = lower(input(text_str, 's'));
    acceptable_inputs = {'yes', 'no'};
    while ~any(strcmp(delete_orig_data, acceptable_inputs))
        warning('Input must be "yes" or "no". Please try again.')
        delete_orig_data = lower(input(text_str, 's'));
    end

    if strcmp(delete_orig_data, 'yes')
        % Delete originally downloaded data
        disp('Deleting data from original download location:');
        [~, delete_message] = rmdir(downloaded_data_loc, 's');
        if isempty(delete_message)
            disp('Finished deleting the original downloaded data!')
        else
            warning(delete_message);
            fprintf('Please manually delete the originally downloaded data located at %s', downloaded_data_loc);
        end
    end

else
    % If the user chose not to use the automated copying functionality, you
    % will need to manually move the data to the correct location in the
    % data folder hierarchy.
    data_path = get_data_path('path_type', 'doppler');
    project_record_path = get_data_path('path_type', 'root');
    disp('-------------')
    disp('You chose to not automatically copy the downloaded data to the correct locations in the data folder.')
    disp('Please move the downloaded Doppler data and Project Record JSON from the CaltechData repository to respective locations below:');
    disp('****************');
    fprintf('\nData path - %s', data_path);
    fprintf('\nJSON path - %s', project_record_path);
    disp('****************\n');

end