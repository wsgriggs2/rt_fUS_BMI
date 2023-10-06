% rigid_body_alignment   Rigid body alignment of two neurovascular images
%
% INPUTS:
%   image1:              2D matrix; Neurovascular map #1
%   image2:              2D matrix; Neurovascular map #2
%   
% OUTPUTS:
%   rigidbody_tform:     (3 x 3) double matrix; Transform matrix
%   image2_tform:        MATLAB format for transform matrix. Either
%                        `affine2d` or `affinetform2d` depending on MATLAB
%                        version.
%   vertshift:           Pixel shift - vertical
%   horzshift:           Pixel shift - horizontal
%   rotation:            Angular rotation
%   image2_transformed   Image2 with transforms applied