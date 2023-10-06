%% Script that will run searchlight analysis and save results


%% User defined settings
split_on_sulcus = true;
remove_sulcus_pixels = true;


%% Load data
data = load_doppler_data('single_or_multiple', 'multiple');

% Spatial filter
data.dop = preprocess_data(data.dop,...
    'spatial_filter', {'disk', 2, 0});

[y_pix, x_pix, n_timepoints] = size(data.dop);


%% Extract train and train_labels
% This also performs the rolling z-score, so no detrend is needed and no
% z-scoring in crossValidate is needed.
[train_data, train_labels] = extract_training_data_and_labels(data);

% Reshape to be (nTrials x features) where features is 1D vector of
% Pixels x time
train_data = reshape(train_data, [], size(train_data, 4))';

% How many classes do we have?
unique_classes = unique(train_labels);
num_classes = length(unique_classes);
switch num_classes
    case 2
        decode_type = '2 tgt';
    case 8
        decode_type = '8 tgt';
    otherwise
        error('This number of classes is not supported.');
        
end

%% Searchlight analysis
% Use the best-time index and ran prediction at single voxel level. Then
% display a map showing prediction ability of individual pixels.
switch decode_type
    case '8 tgt'
        % Takes a while to generate.
        num_replicates = 1000; % For the angular error
    case '2 tgt'
        num_replicates = 1000; % For the angular error
end
validation_string = 'kFold';


% choose K for K-fold validation
K = 3;

radii_to_test = [2];

[cols_in_img, rows_in_img] = meshgrid(1:x_pix, 1:y_pix);


% For speed purposes, load and pass the null distribution.
% To make smaller, we can shorten to length of train_labels. This
% will still be longer than necessary, but sets an upper bound on
% possible # of trials.
n_trials = length(train_labels);
null_mean_angular_error = load_null_mean_angular_error(num_classes, num_replicates, n_trials);

if split_on_sulcus || remove_sulcus_pixels
    if size(data.session_run_list, 1) > 1
        warning('Sulcus and brain identification: Assuming that the anatomy will be the same throughout all concatenated data.');
    end
    sulcus_vertices = get_sulcus_info(data.session_run_list(1, 1), data.session_run_list(1, 2), data.neurovascular_map, data.UF, ...
        'verbose', true);
    
    sulcus_pixels = poly2mask(sulcus_vertices(:, 1), sulcus_vertices(:, 2), y_pix, x_pix);
end

switch decode_type
    case '8 tgt'
        classifier_string = 'PCA+LDA';
        m = 7;
        
        
        lookup_table_for_horz = [1 2 3 1 2 3 1 2 3];
        lookup_table_for_vert = [1 1 1 2 2 2 3 3 3];
        train_labels_horz = lookup_table_for_horz(train_labels);
        train_labels_vert = lookup_table_for_vert(train_labels);
        train_labels_combined = [train_labels_horz' train_labels_vert'];
        
        [sl_percentCorrect_horz, sl_percentCorrect_vert, ...
            sl_percentCorrect_combined, sl_pvalue_combined,...
            sl_angularError_combined, sl_angularError_pvalue_combined] = deal(NaN(y_pix, x_pix, length(radii_to_test)));
        
        [sl_confusion_horz, sl_confusion_vert, sl_confusion_combined, sl_result_horz, sl_result_vert, sl_result_combined] = deal(cell(y_pix, x_pix, length(radii_to_test)));
        
        p = gcp; % Initialize the parallel pool.
        fprintf('Searchlight analysis started at: %s\n', datestr(now()));
        for k = 1:length(radii_to_test)
            radius = radii_to_test(k);
            fprintf('Searchlight analysis for radius %0.2f ... ', radius);
            
            ppm = ParforProgressbar(y_pix * x_pix);
            % Requires Parallel Computing toolbox. If not available, then
            % switch to `for` loop.
            parfor i = 1:x_pix
                for j = 1:y_pix
                    if (~split_on_sulcus && ~remove_sulcus_pixels) || ~sulcus_pixels(j, i)
                        % prepare the training data for this time window
                        ROI_pixels = (rows_in_img - j).^2 ...
                            + (cols_in_img - i).^2 <= radius.^2;
                        
                        
                        
                        % If sulcus splits circle pixels, then only keep half with the
                        % center pixel.
                        if split_on_sulcus
                            ROI_sulcus_overlap = ROI_pixels & sulcus_pixels;
                            if nnz(ROI_sulcus_overlap)
                                ROI_sulcus_exclude = ROI_pixels & ~sulcus_pixels;
                                binary_comps = bwlabel(ROI_sulcus_exclude);
                                center_pixel_comp = binary_comps(j, i);
                                ROI_pixels = binary_comps == center_pixel_comp;
                            end
                        elseif remove_sulcus_pixels
                            ROI_pixels = ROI_pixels & ~sulcus_pixels;
                        end
                        
                        mask_ind = false(y_pix, x_pix);
                        mask_ind(ROI_pixels) = true;
                        
                        mask_ind_1D = reshape(mask_ind, 1, []);
                        mask_ind_1D = repmat(mask_ind_1D, 1, 3);
                        
                        train_data_sl = train_data(:,mask_ind_1D);
                        
                        
                        % cross validated performance is calculated here
                        [cp_horz, p_horz, cp_vert, p_vert, cp_combined, p_combined, custom_counting_matrix] = cross_validate_multicoder(train_data_sl, train_labels_combined, ...
                            'classificationMethod', classifier_string, ...
                            'validationMethod', validation_string, ...
                            'K', K);
                        
                        % Angular error calculation
                        % Reconstruct a record for final performance. This
                        % is not the actual order of predictions since that
                        % is not easily accessible. This is correct number
                        % of each prediction though.
                        [actual_class, predicted_class] = deal([]);
                        for true_class_num = 1:9
                            for predicted_class_num = 1:9
                                actual_class = [actual_class; true_class_num*ones(custom_counting_matrix(predicted_class_num,true_class_num), 1)];
                                predicted_class = [predicted_class; predicted_class_num*ones(custom_counting_matrix(predicted_class_num, true_class_num), 1)];
                            end
                        end
                        [empiric_mean_angular_error, ~, pvalue_empiric] = calculate_angular_error(predicted_class, actual_class, ...
                            'num_replicates', num_replicates, ...
                            'null_mean_angular_error', null_mean_angular_error);
                        
                        % grabbing useful values
                        sl_percentCorrect_horz(j,i, k) = cp_horz.CorrectRate*100;
                        sl_confusion_horz{j,i, k} = cp_horz.CountingMatrix;
                        
                        sl_percentCorrect_vert(j,i, k) = cp_vert.CorrectRate*100;
                        sl_confusion_vert{j,i, k} = cp_vert.CountingMatrix;
                        
                        sl_percentCorrect_combined(j,i, k) = cp_combined.CorrectRate*100;
                        sl_pvalue_combined(j, i, k) = p_combined;
                        sl_confusion_combined{j,i, k} = custom_counting_matrix;
                        sl_angularError_combined(j, i, k) = empiric_mean_angular_error;
                        sl_angularError_pvalue_combined(j, i, k) = pvalue_empiric;
                        
                        % save the data to results
                        sl_result_horz{j,i, k}.percentCorrect =sl_percentCorrect_horz(j,i,k);
                        sl_result_horz{j,i, k}.confusion = sl_confusion_horz{j,i,k};
                        sl_result_horz{j,i, k}.p = p_horz;
                        
                        sl_result_vert{j,i, k}.percentCorrect = sl_percentCorrect_vert(j,i,k);
                        sl_result_vert{j,i, k}.confusion = sl_confusion_vert{j,i,k};
                        sl_result_vert{j,i, k}.p = p_vert;
                        
                        sl_result_combined{j,i, k}.percentCorrect = sl_percentCorrect_combined(j,i,k);
                        sl_result_combined{j,i, k}.confusion = sl_confusion_combined{j,i,k};
                        sl_result_combined{j,i, k}.p = p_combined;
                        ppm.increment();
                    end
                end
            end
            delete(ppm);
            fprintf('Completed %s\n', datestr(now()));
        end
        
    case '2 tgt'
        
        classifier_string = 'CPCA+LDA';
        m = 1;
        
        [sl_percentCorrect_combined, sl_pvalue_combined] = deal(NaN(y_pix, x_pix, length(radii_to_test)));
        
        [sl_confusion_combined, sl_result_combined] = deal(cell(y_pix, x_pix, length(radii_to_test)));
        
        if size(train_labels, 1) == 1 && size(train_labels, 2) > 1
            train_labels = train_labels';
        end
        
        p = gcp; % Initialize the parallel pool.
        fprintf('Searchlight analysis started at: %s\n', datestr(now()));
        for k = 1:length(radii_to_test)
            radius = radii_to_test(k);
            fprintf('Searchlight analysis for radius %0.2f ... ', radius);
            
            ppm = ParforProgressbar(y_pix * x_pix);
            % Requires Parallel Computing toolbox. If not available, then
            % switch to `for` loop.
            parfor i = 1:x_pix
                for j = 1:y_pix
                    if (~split_on_sulcus && ~remove_sulcus_pixels) || ~sulcus_pixels(j, i)
                        % prepare the training data for this time window
                        ROI_pixels = (rows_in_img - j).^2 ...
                            + (cols_in_img - i).^2 <= radius.^2;
                        
                        
                        
                        % If sulcus splits circle pixels, then only keep half with the
                        % center pixel.
                        if split_on_sulcus
                            ROI_sulcus_overlap = ROI_pixels & sulcus_pixels;
                            if nnz(ROI_sulcus_overlap)
                                ROI_sulcus_exclude = ROI_pixels & ~sulcus_pixels;
                                binary_comps = bwlabel(ROI_sulcus_exclude);
                                center_pixel_comp = binary_comps(j, i);
                                ROI_pixels = binary_comps == center_pixel_comp;
                            end
                        elseif remove_sulcus_pixels
                            ROI_pixels = ROI_pixels & ~sulcus_pixels;
                        end
                        
                        mask_ind = false(y_pix, x_pix);
                        mask_ind(ROI_pixels) = true;
                        
                        mask_ind_1D = reshape(mask_ind, 1, []);
                        mask_ind_1D = repmat(mask_ind_1D, 1, 3);
                        
                        train_data_sl = train_data(:,mask_ind_1D);
                        
                        % cross validated performance is calculated here
                        [cp_combined, p_combined] = cross_validate(train_data_sl, train_labels, ...
                            'classificationMethod', classifier_string,  ...
                            'validationMethod', validation_string,...
                            'K', K);
                        
                        % Angular error calculation
                        % Reconstruct a record for final performance. This
                        % is not the actual order of predictions since that
                        % is not easily accessible. This is correct number
                        % of each prediction though.
                        [actual_class, predicted_class] = deal([]);
                        counting_matrix = cp_combined.CountingMatrix;
                        for true_class_num = 1:2
                            for predicted_class_num = 1:2
                                actual_class = [actual_class; true_class_num*ones(counting_matrix(predicted_class_num,true_class_num), 1)];
                                predicted_class = [predicted_class; predicted_class_num*ones(counting_matrix(predicted_class_num, true_class_num), 1)];
                            end
                        end
                        
                        % Null_mean_angular_error is loaded once before the
                        % searchlight is run for speed reasons.
                        [empiric_mean_angular_error, ~, pvalue_empiric] = calculate_angular_error(predicted_class, actual_class, ...
                            'num_replicates', num_replicates, ...
                            'null_mean_angular_error', null_mean_angular_error);
                        
                        
                        % grabbing useful values
                        
                        sl_percentCorrect_combined(j,i, k) = cp_combined.CorrectRate*100;
                        sl_pvalue_combined(j, i, k) = p_combined;
                        sl_confusion_combined{j,i, k} = cp_combined.CountingMatrix;
                        sl_angularError_combined(j, i, k) = empiric_mean_angular_error;
                        sl_angularError_pvalue_combined(j, i, k) = pvalue_empiric;
                        
                        % save the data to result
                        
                        sl_result_combined{j,i, k}.percentCorrect = sl_percentCorrect_combined(j,i,k);
                        sl_result_combined{j,i, k}.confusion = sl_confusion_combined{j,i,k};
                        sl_result_combined{j,i, k}.p = p_combined;
                        ppm.increment();
                    end
                end
            end
            delete(ppm);
            fprintf('Completed %s\n', datestr(now()));
        end
end % End for switch statement for searchlighy analysis for different classes (2 vs 8)

%% Save data to file (if desired)
if size(data.session_run_list, 1) > 1
    suggested_fname = sprintf('searchlight_results_multipleSessions_gen%s.mat',...
        datestr(now,'yyyymmdd'));
else
    suggested_fname = sprintf('searchlight_results_S%dR%d_gen%s.mat',...
        data.session_run_list(:, 1), data.session_run_list(:, 2), ...
        datestr(now,'yyyymmdd'));
end

title_string = 'Specify where to save searchlight results';
disp(title_string);
[filename, filepath] = uiputfile(fullfile(get_data_path('path_type', 'output'), suggested_fname), title_string);


neurovascular_map = data.neurovascular_map;
UF = data.UF;
session_run_list = data.session_run_list;
if any([filename filepath])
    save(fullfile(filepath, filename), 'neurovascular_map', 'UF', ...
        'sl_pvalue_combined', 'sl_angularError_pvalue_combined', ...
        'sl_angularError_combined', 'decode_type', 'sl_percentCorrect_combined', ...
        'radii_to_test', 'session_run_list');
end

%% Visualize searchlight analysis

display_accuracy_or_angular_error = 'angular_error'; % 'angular_error' or 'accuracy'
display_type = 'bottom_voxels'; % 'bottom_voxels' or 'top_voxels' Or 'pvalue_threshold'

percent_of_voxels_to_keep = .1; % In decimal form, so 0.1 = 10%
pvalue_threshold = 1e-2;
pixelsize = 0.1;
X_img_mm = pixelsize/2 + (0:x_pix-1)*pixelsize; % So that the center of the image is at true halfway point, not shifted.
Z_img_mm = pixelsize/2 + (0:y_pix-1)*pixelsize + data.UF.Depth(1);


figure;
tld = tiledlayout('flow');
circle_center = [ 2 max(Z_img_mm)-2];

% Apply FDR correction for the p-values
pvalues_1D = reshape(sl_pvalue_combined,[],1);
nan_ind = isnan(pvalues_1D);
[Q] = mafdr(pvalues_1D(~nan_ind), 'BHFDR', true);
pvalues_corrected = NaN(size(pvalues_1D));
pvalues_corrected(~nan_ind) = Q;
matrix_size = size(sl_pvalue_combined);
sl_pvalue_combined_corrected = reshape(pvalues_corrected, matrix_size);

% Apply FDR correction for the angular error p-values
pvalues_1D = reshape(sl_angularError_pvalue_combined,[],1);
nan_ind = isnan(pvalues_1D);
[Q] = mafdr(pvalues_1D(~nan_ind), 'BHFDR', true);
pvalues_corrected = NaN(size(pvalues_1D));
pvalues_corrected(~nan_ind) = Q;
matrix_size = size(sl_angularError_pvalue_combined);
sl_angularError_pvalue_combined_corrected = reshape(pvalues_corrected, matrix_size);

%Ignore voxels that are greater than pi/2 (chance level) for angular error
sl_angularError_combined_corrected = sl_angularError_combined;
% above_pi_half_ind = sl_angularError_combined>=pi/2;
% sl_angularError_combined_corrected(above_pi_half_ind) = NaN;
% sl_angularError_pvalue_combined_corrected(above_pi_half_ind) = NaN;

switch decode_type
    case '2 tgt'
        colormap_to_use = inferno;
    case '8 tgt'
        colormap_to_use = viridis;
end

switch display_accuracy_or_angular_error
    case 'angular_error'
        pvalue_to_use = sl_angularError_pvalue_combined_corrected;
        performance_metric_to_use = sl_angularError_combined_corrected * 180/pi;
        colorbar_title = 'Angular error (deg)';
        unit_measure = ' deg';
        colorscale = flipud(colormap_to_use);
        
        colorbar_max = 90;
        colorbar_min = 30;
    case 'accuracy'
        pvalue_to_use = sl_pvalue_combined_corrected;
        performance_metric_to_use = sl_percentCorrect_combined;
        colorbar_title = 'Performance (% correct)';
        unit_measure = '%';
        colorscale = colormap_to_use;
        
        colorbar_max = 100;
        colorbar_min = 0;
end


colorbar_limits = [0 max(performance_metric_to_use,[], 'all')]; % Same color range for all searchlight analysis plots, scaled to max performance.


for k = 1:length(radii_to_test)
    
    % Plot searchlight analysis results
    
    
    % Overlay performance on anatomy
    switch display_type
        case 'pvalue_threshold'
            pvalue_mask = pvalue_to_use < pvalue_threshold;
            
            performance_masked = squeeze(performance_metric_to_use(:, :, k));
            performance_masked(~pvalue_mask) = NaN;
            title_string = sprintf('thresholded at p<=%d', pvalue_threshold);
        case 'top_voxels'
            % Find top X% of voxels
            performance = performance_metric_to_use(:, :, k);
            pvalue = pvalue_to_use(:,:, k);
            top_quantile_threshold(k) = quantile(performance(:), 1-percent_of_voxels_to_keep);
            performance_mask = performance >= top_quantile_threshold(k);
            performance_masked = performance;
            performance_masked(~performance_mask) = NaN;
            title_string = sprintf('keeping top %d%% of voxels', percent_of_voxels_to_keep*100);
            
            
            
            % Report the p-value associated with the respective quantile thresholds
            pvalue_of_quantile_threshold = max(pvalue(performance <= top_quantile_threshold(k)),[], 'all');
            fprintf('Radius - %d p-value associated with the quantile threshold of %0.2f%s correct is %d\n', radii_to_test(k), top_quantile_threshold(k), unit_measure, pvalue_of_quantile_threshold);
        case 'bottom_voxels'
            % Find top X% of voxels
            performance = performance_metric_to_use(:, :, k);
            pvalue = pvalue_to_use(:,:, k);
            bottom_quantile_threshold(k) = quantile(performance(:), percent_of_voxels_to_keep);
            performance_mask = performance <= bottom_quantile_threshold(k);
            performance_masked = performance;
            performance_masked(~performance_mask) = NaN;
            title_string = sprintf('keeping top %d%% of voxels', percent_of_voxels_to_keep*100);
            
            % Report the p-value associated with the respective quantile thresholds
            pvalue_of_quantile_threshold = max(pvalue(performance <= bottom_quantile_threshold(k)),[], 'all');
            fprintf('Radius - %d p-value associated with the quantile threshold of %0.2f%s correct is %d\n', radii_to_test(k), bottom_quantile_threshold(k), unit_measure, pvalue_of_quantile_threshold);
    end
    
    % Update colorbar limits
    colorbar_max = max([colorbar_max; performance_masked(:)], [], 'all');
    colorbar_min = min([colorbar_min; performance_masked(:)], [], 'all');
    
    performance_masked_acrossRadii{k} = performance_masked;
end
for k=1:length(radii_to_test)
    radius = radii_to_test(k);
    nexttile(k);
    
    background_image = data.neurovascular_map;
    
    cmap = plotDuplexImage(X_img_mm, Z_img_mm, performance_masked_acrossRadii{k}, background_image,...
        'colormap2use', colorscale, 'nonlinear_bg',2, ...
        'AutoColorBarLimits', [colorbar_min colorbar_max], ...
        'showColorbar', true, ...
        'ColorbarTitle', colorbar_title);
    if size(data.session_run_list, 1)>1
        session_string = sprintf('%d ', data.session_run_list(:, 1)');
        title(sprintf('(radius = %0.2f) \nMultiple sessions - [%s]', radius, session_string));
    else
        title(sprintf('(radius = %0.2f) \nS%dR%d', radius, data.session_run_list(:, 1), data.session_run_list(:, 2)));
    end
    hold on;
    circle_handle = viscircles(circle_center,radius*pixelsize, 'EnhanceVisibility', false, ...
        'Color', 'w');
end

title(tld, sprintf('Searchlight analysis %s', title_string));





