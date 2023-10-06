function [cp, p] = cross_validate(data, labels, varargin)
% cross_validate  Perform cross-validated decoding using model of your
% choice.
%
% INPUTS:
%   data:                       2D double; [nImages x (xPixels*yPixels)]; 
%                               training data
%   labels:                     (1 x nTrials) double; Label associated with
%                               each sample
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
%   cp:                         'class performance', created using classperf
%                               note that, among other things, cp contains:
%                               confusion matrix: 'countingMatrix', e.g. 
%                               s[1 0 ; 0 1]
%                               correct rate, e.g. 0.95
%   p:                          p value (based on binomial test), e.g. 0.05
%
% See also cross_validate_multicoder

%% handling varargin
inp = inputParser;
inp.addOptional('verbose',false,@islogical)
inp.addOptional('validationMethod','kFold')
inp.addOptional('K',10);
inp.addOptional('classificationMethod','CPCA+LDA')
inp.addParameter('m', 1);
inp.addParameter('variance_to_keep', 95);
inp.addOptional('trial_ind', 1:length(labels));
inp.parse(varargin{:});
sets = inp.Results;
m = sets.m;

%% get useful vars
N = size(data,1);           % number of trials
cp = classperf(labels);     % initializing class performance var

%% k-fold & leave one out cross validation
    
    % leave one out is just a special case of K-Fold wehre K=N;
    if strcmp(sets.validationMethod,'leaveOneOut')
        sets.K = N;
    end
    
    % cross validate here
    indices = crossvalind('Kfold', N, sets.K);
    % for each k-fold
    if sets.verbose
        fprintf('%d folds: ', sets.K);
    end
    for i = 1:sets.K
        % create indices for training set & test set
        test = (indices == i);
        train = ~test;
        
        %If nans, then define as 0. NaN values will break downstream
        %functions.
        data(isnan(data)) = 0;
        
        % Train model
        model = train_classifier(data(train, :), labels(train, :), ...
            'method', sets.classificationMethod, ...
            'variance_to_keep', sets.variance_to_keep);
        
        % Test model on held-out data
        class = make_prediction(data(test, :), model);
        
        % Assess performance
        classperf(cp, class, test);
        
        % Update command line output (if requested)
        if sets.verbose
            fprintf('|');
        end
    end
    if sets.verbose
        fprintf('\n');
    end


%% calculate classification accuracy measures
percentCorrect = cp.correctRate*100;
nCorrect = sum(diag(cp.CountingMatrix));
nCounted = sum(cp.CountingMatrix(:));
chance = 1/length(unique(labels));
p = binomialTest(nCorrect, nCounted, chance, 'one');

%% display measures if verbose is on
if sets.verbose
    % classification accuracy (%)
    fprintf('\nClassification Accuracy: \n%i / %i trials correctly classified (%2.2f%% correct)\t',...
        nCorrect, nCounted, percentCorrect)
end

end
