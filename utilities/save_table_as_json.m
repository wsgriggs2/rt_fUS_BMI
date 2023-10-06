function success = save_table_as_json(table, filename)

fid = fopen(filename, 'w');
if verLessThan('matlab', '9.1')
    error('This version of MATLAB is too old. Please update to at least MATLAB 2016b to have JSON encoding.');  
elseif verLessThan('matlab', '9.10')
    encoded_json = jsonencode(table);
    warning('Saved JSON file, but please upgrade to MATLAB 2021a for new JSON capabilities');
else
    % Added pretty print capability in 2021a (ver 9.10)
    encoded_json = jsonencode(table, 'PrettyPrint', true);
end
success = fwrite(fid, encoded_json);
fclose(fid);

