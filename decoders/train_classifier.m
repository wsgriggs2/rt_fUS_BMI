% train_classifier Train the specified classifier
%
% INPUTS:
%   trainData:              (nSamples x yPix*xPix*nTimepoints); 2D array
%   trainLabels:            (1 x nTrials) double; Label associated with
%                           each sample
%   varargin
%       'method':           string; Accepted values are:
%                           'CPCA+LDA',
%                           'PCA+LDA',
%       'variance_to_keep`: double between (0 100); What percent of the
%                           variance do you want to keep?
%       'm':                The dimensionality of the CPCA subspaces. The
%                           max is (nClasses-1).
%
% OUTPUTS:
%   model:                  The trained model and relevant information
%                           about the model.
%
% See also make_prediction
