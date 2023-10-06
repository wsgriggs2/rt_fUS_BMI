function [cp_horz, p_horz, cp_vert, p_vert, cp_combined, p_combined, countingmatrix_custom] = cross_validate_multicoder(data, labels, varargin)
% cross_validate_multicoder  Perform cross-validated decoding using model
%                            of your choice. Uses `multicoder` architecture
%                            where we separately decode horizontal and
%                            vertical, then combine predictions.
%
% INPUTS:
%   data:                       2D double; [nImages x (xPixels*yPixels)];
%                               training data
%   labels:                     (2 x nTrials) double; Label associated with
%                               each sample. x, y columns
%                               1 = Negative (down or left)
%                               2 = At center
%                               3 = Positive (up or right)
%   varargin:
%       verbose:                bool; print progress to command line?
%       validationMethod:       string; Method for cross-validation.
%                               Options are: 'leaveOneOut' or 'kFold'
%       K:                      scalar double; If using 'kFold method',
%                               then how many folds?
%       classificationMethod:   String; Accepted values are: 'CPCA+LDA' or
%                               'PCA+LDA'
%       m:                      Scalar positive value; determines dimension
%                               of each cPCA subspace. Max is (# of classes
%                               - 1).
%       variance_to_keep:       scalar double; Percent of variance to keep;
%                               Out of 100. For the PCA dimensionality
%                               reduction methods. Not currently used for
%                               cPCA.
%       trial_ind:
%
% OUTPUTS:
%   countingmatrix_custom:      Confusion matrix for the combined
%                               predictions.s
%
%   Separately for horizontal, vertical, and combined
%   cp:                         'class performance', created using classperf
%                               note that, among other things, cp contains:
%                               confusion matrix: 'countingMatrix', e.g.
%                               s[1 0 ; 0 1]
%                               correct rate, e.g. 0.95
%   p:                          p value (based on binomial test), e.g. 0.05
%
% See also cross_validate

%

%% handling varargin
p = inputParser;
p.addOptional('verbose',false,@islogical)
p.addOptional('validationMethod','kFold')
p.addOptional('K',10);
p.addOptional('classificationMethod','CPCA+LDA')
p.addParameter('m', 1)
p.addParameter('variance_to_keep', 95);
p.parse(varargin{:});
sets = p.Results;

%% get useful vars
display_progress = sets.verbose;
N = size(data,1);           % number of trials

cp_horz = classperf(labels(~isnan(labels(:,1)),1));     % initializing class performance var
cp_vert = classperf(labels(~isnan(labels(:,2)),2));     % initializing class performance var

unique_horz_classes = unique(labels(~isnan(labels(:, 1))));
unique_vert_classes = unique(labels(~isnan(labels(:, 2))));

possible_combined_labels = length(unique_horz_classes) * length(unique_vert_classes);

countingmatrix_custom = zeros(possible_combined_labels);

labels_combined = labels(:,1) + length(unique(labels(~isnan(labels(:,1)),1))) * (labels(:,2)-1);

if all(~isnan(labels_combined))
    cp_combined = classperf(labels_combined);
end

%% k-fold & leave one out cross validation


% leave one out is just a special case of K-Fold wehre K=N;
if strcmp(sets.validationMethod,'leaveOneOut')
    sets.K = N;
    display_progress  = true;
end

% cross validate here
indices = crossvalind('Kfold', N, sets.K);
% for each k-fold
for i = 1:sets.K
    % create indices for training set & test set
    test = (indices == i);
    train = ~test;
    nan_ind = isnan(labels);
    
    % Accounting for NaN labels in case we are not using all possible
    % classes of data.
    horz_train = train & ~nan_ind(:,1);
    vert_train = train & ~nan_ind(:,2);
    horz_test = test & ~nan_ind(:,1);
    vert_test = test & ~nan_ind(:,2);
    
    % Fit and apply z-scoring to train data, then apply to test.
    [data(train, :), mu, sigma] = zscore(data(train, :));
    data(test, :) = (data(test, :) - mu) ./ sigma;
    
    %If nans, then define as 0. NaN values will break downstream
    %functions.
    data(isnan(data)) = 0;
    
    % Train models
    model_horz = train_classifier(data(horz_train, :), ...
        labels(horz_train,1), ...
        'method', sets.classificationMethod, ...
        'variance_to_keep', sets.variance_to_keep);
    model_vert = train_classifier(data(vert_train, :), ...
        labels(vert_train,2), ...
        'method', sets.classificationMethod, ...
        'variance_to_keep', sets.variance_to_keep);
    
    % Test models using held-out data
    class_horz = make_prediction(data(horz_test, :), model_horz);
    class_vert = make_prediction(data(vert_test, :), model_vert);
    
    if all(~isnan(labels_combined))
        class_combined = class_horz + length(unique(labels(~isnan(labels(:,1)),1))) * (class_vert-1);
    end
    
    classperf(cp_horz, class_horz, horz_test(~isnan(labels(:,1))));
    classperf(cp_vert, class_vert, vert_test(~isnan(labels(:,2))));
    if all(~isnan(labels_combined))
        classperf(cp_combined, class_combined, test);
        
        countingmatrix_custom = update_counting_matrix(countingmatrix_custom, class_combined, labels_combined(test));
    end
    if display_progress
        if i == 1
            f = waitbar(0, 'Leave-one-out analysis - Please wait...');
        else
            waitbar(i/sets.K, f);
        end
    end
end
if display_progress
    close(f);
end


%% calculate classification accuracy measures
% May be slight mismatch between percentCorrect and nCorrect/nCounted due
% to `inconclusive` samples. MATLAB ignores them. We keep them in nCounted.
percentCorrect_horz = cp_horz.correctRate*100;
nCorrect_horz = sum(diag(cp_horz.CountingMatrix));
nCounted_horz = sum(cp_horz.CountingMatrix(:));
chance_horz = 1/length(cp_horz.ClassLabels);
p_horz = binomialTest(nCorrect_horz, nCounted_horz, chance_horz, 'one');

percentCorrect_vert = cp_vert.correctRate*100;
nCorrect_vert = sum(diag(cp_vert.CountingMatrix));
nCounted_vert = sum(cp_vert.CountingMatrix(:));
chance_vert = 1/length(cp_vert.ClassLabels);
p_vert = binomialTest(nCorrect_vert, nCounted_vert, chance_vert, 'one');

percentCorrect_combined = cp_combined.correctRate*100;
nCorrect_combined = sum(diag(cp_combined.CountingMatrix));
nCounted_combined = sum(cp_combined.CountingMatrix(:));
chance_combined = 1/length(cp_combined.ClassLabels);
p_combined = binomialTest(nCorrect_combined, nCounted_combined, chance_combined, 'one');

%% display measures if verbose is on
if sets.verbose
    % classification accuracy (%)
    fprintf('\nHorizontal Classification Accuracy: \n%i / %i trials correctly classified (%2.2f%% correct)\t',...
        nCorrect_horz, nCounted_horz, percentCorrect_horz)
    
    % classification accuracy (%)
    fprintf('\nVertical Classification Accuracy: \n%i / %i trials correctly classified (%2.2f%% correct)\t',...
        nCorrect_vert, nCounted_vert, percentCorrect_vert)
    
    % classification accuracy (%)
    fprintf('\nCombined Classification Accuracy: \n%i / %i trials correctly classified (%2.2f%% correct)\t',...
        nCorrect_combined, nCounted_combined, percentCorrect_combined)
end

end
