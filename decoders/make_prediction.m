% make_prediction  Take the trained model and predict the intended class.
% This takes in the test data and previously trained model and returns the
% predicted class.
%
% INPUTS:
%   testData:   (1 x yPix*xPix*nTimepoints); A single test data point
%   model:      struct; Previously trained model; can be a variety of types,
%               including `PCA+LDA` and 'CPCA+LDA'.
%      
% OUTPUTS:
%   class:      (int); Predicted class
%
% See also train_classifier

