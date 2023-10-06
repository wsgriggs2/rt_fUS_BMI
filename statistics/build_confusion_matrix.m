function confusion_mat = build_confusion_matrix(n_classes, actual, predicted)
% build_confusion_matrix  Build a confusion matrix
%
% INPUTS:
%   n_classes:      scalar double; How many different classes
%   actual:         (n x 1); Actual labels
%   predicted:      (n x 1); Predited labels
%      
% OUTPUTS:
%   confusion_mat:  (n_classes x n_classes) double; Matrix of counts for
%                   each predict-actual class pair


confusion_mat = zeros(n_classes, n_classes);
for ac = 1:n_classes
    for pr = 1:n_classes
        confusion_mat(ac, pr) = nnz(actual == ac & predicted == pr);
    end
end
end