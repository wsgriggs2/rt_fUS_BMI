%% Script to visualize searchlight analysis results

%% Load data
suggested_fname = 'searchlight_results_S*R*_gen*.mat';
title_string = 'Select searchlight results to visualize.';
disp(title_string);
[filename, filepath] = uigetfile(fullfile(get_data_path('path_type', 'output'), ...
    suggested_fname), title_string);

data = load(fullfile(filepath, filename), 'neurovascular_map', 'UF', ...
        'sl_pvalue_combined', 'sl_angularError_pvalue_combined', ...
        'sl_angularError_combined', 'decode_type', 'sl_percentCorrect_combined', ...
        'radii_to_test', 'session_run_list');

%% Generate some basic variables needed
[y_pix, x_pix] = size(data.neurovascular_map);

%%

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
pvalues_1D = reshape(data.sl_pvalue_combined,[],1);
nan_ind = isnan(pvalues_1D);
[Q] = mafdr(pvalues_1D(~nan_ind), 'BHFDR', true);
pvalues_corrected = NaN(size(pvalues_1D));
pvalues_corrected(~nan_ind) = Q;
matrix_size = size(data.sl_pvalue_combined);
sl_pvalue_combined_corrected = reshape(pvalues_corrected, matrix_size);

% Apply FDR correction for the angular error p-values
pvalues_1D = reshape(data.sl_angularError_pvalue_combined,[],1);
nan_ind = isnan(pvalues_1D);
[Q] = mafdr(pvalues_1D(~nan_ind), 'BHFDR', true);
pvalues_corrected = NaN(size(pvalues_1D));
pvalues_corrected(~nan_ind) = Q;
matrix_size = size(data.sl_angularError_pvalue_combined);
sl_angularError_pvalue_combined_corrected = reshape(pvalues_corrected, matrix_size);



switch data.decode_type
    case '2 tgt'
        colormap_to_use = inferno;
    case '8 tgt'
        colormap_to_use = viridis;
end

switch display_accuracy_or_angular_error
    case 'angular_error'
        pvalue_to_use = sl_angularError_pvalue_combined_corrected;
        performance_metric_to_use = data.sl_angularError_combined * 180/pi;
        colorbar_title = 'Angular error (deg)';
        unit_measure = ' deg';
        colorscale = flipud(colormap_to_use);
        
        colorbar_max = 90;
        colorbar_min = 30;
    case 'accuracy'
        pvalue_to_use = sl_pvalue_combined_corrected;
        performance_metric_to_use = data.sl_percentCorrect_combined;
        colorbar_title = 'Performance (% correct)';
        unit_measure = '%';
        colorscale = colormap_to_use;
        
        colorbar_max = 100;
        colorbar_min = 0;
end


colorbar_limits = [0 max(performance_metric_to_use,[], 'all')]; % Same color range for all searchlight analysis plots, scaled to max performance.


for k = 1:length(data.radii_to_test)
    
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
            fprintf('Radius - %d p-value associated with the quantile threshold of %0.2f%s correct is %d\n', data.radii_to_test(k), top_quantile_threshold(k), unit_measure, pvalue_of_quantile_threshold);
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
            fprintf('Radius - %d p-value associated with the quantile threshold of %0.2f%s correct is %d\n', data.radii_to_test(k), bottom_quantile_threshold(k), unit_measure, pvalue_of_quantile_threshold);
    end
    
    % Update colorbar limits
    colorbar_max = max([colorbar_max; performance_masked(:)], [], 'all');
    colorbar_min = min([colorbar_min; performance_masked(:)], [], 'all');
    
    performance_masked_acrossRadii{k} = performance_masked;
end
for k=1:length(data.radii_to_test)
    radius = data.radii_to_test(k);
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