function plt = create_confusion_matrix(counting_matrix, varargin)
% create_confusion_matrix  Create a confusion matrix using imagesc
% X-axis is predicted value
% Y-axis is true value
%
% INPUTS:
% counting_matrix:          (n x m) 2D matrix; Each row reflects the predicted
%                           guesses for that true label. Each column represents
%                           the true labels for a specific prediction.
%
%   varargin: 
%       label_strings:      cell of strings; If used, Must have one label for 
%                           each predicted/true label, i.e. must have max(n,m)
%                           labels.
%       row_normalize:      bool; Do you want each row to sum to 1?
%       title:              string; Title for the plot
%       FontColor:          string or 1x3 RGB array;
%       superimpose_text:   string; either 'count' or 'percentage'; if empty,
%                           nothing will be superimposed
%       colorbar_limits:    1x2 array with the percentage limits for the colorbar
%       colormap:           Colormap that you want to use for imagesc. Default
%                           is inferno.
%      
% OUTPUTS:
%   plt:                Handle to the imagesc graphic
%
% See also imagesc


%% Variable input parsing
p = inputParser;
p.addOptional('label_strings', []);
p.addOptional('row_normalize', true);
p.addOptional('title', []);
p.addOptional('FontColor', 'k');
p.addOptional('superimpose_text', '');
p.addOptional('colorbar_limits', [0 100]);
p.addOptional('colormap', inferno);
p.parse(varargin{:});
inputs = p.Results;

%% Identify how many true vs predicted labels
n_true_classes = size(counting_matrix, 1);
n_pred_classes = size(counting_matrix, 2);

%% row normalize
if inputs.row_normalize
    counting_matrix_plot = counting_matrix;
    counting_matrix_plot = counting_matrix_plot./sum(counting_matrix_plot, 2);
    counting_matrix_plot(isnan(counting_matrix_plot)) = 0;
else
    counting_matrix_plot = counting_matrix;
end


%% Create confusion matrix
plt = imagesc(counting_matrix_plot*100);

%% Add items to the imagesc graphic
% Make x-ticks and y-ticks for each predicted and true class
yticks(1:n_true_classes);
xticks(1:n_pred_classes);

% Add labels for each x- and y-tick if requested
if ~isempty(inputs.label_strings)
    yticklabels(inputs.label_strings(1:n_true_classes));
    xticklabels(inputs.label_strings(1:n_pred_classes));
end

% Display colorbar
clb = colorbar;
clb.Label.String = 'Percent (row-normalized)';
clb.Color = inputs.FontColor;

% Set colormap
colormap(inputs.colormap);

% Specify colorbar limits
caxis(inputs.colorbar_limits);

% Label axes
xlabel('Predicted class');
ylabel('True class');

% Add title
title(inputs.title, 'Color', inputs.FontColor);

% Set text color
set(gca, 'XColor', inputs.FontColor);
set(gca, 'YColor', inputs.FontColor);
text_color = inputs.FontColor;

%% Superimpose text on each square of confusion matrix (if requested)
if inputs.superimpose_text
    for i = 1:n_true_classes
        for j = 1:n_pred_classes
            count = counting_matrix(i, j); 
            row_total = sum(counting_matrix(i, :));
            
            % Generate label
            % Can superimpose absolute count or percentage
            switch inputs.superimpose_text
                case 'count'
                    label = strcat(num2str(count), '/', num2str(row_total));
                case 'percentage'
                    label = sprintf('%0.1f%%', 100*counting_matrix_plot(i, j));
                    if counting_matrix_plot(i, j) < .5
                        text_color = 'w';
                    else
                        text_color = 'k';
                    end
            end
            
            % Apply label
            text(j,i,label, ...
                'HorizontalAlignment', 'center', ...
                'Color', text_color)
        end
    end   
end

