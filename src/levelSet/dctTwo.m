function X = dctTwo(X, use_cuda)
% DCTTwo - Perform Discrete Cosine Transform along two dimensions
%
% Syntax: X = dctTwo(X)
%
% Inputs:
%   X - Input data to be transformed
%
% Outputs:
%   X - Transformed data
%
% The function performs the discrete cosine transform along two dimensions.
% The function also casts the input to the GPU to ensure computation takes 
% place on the GPU.

if nargin < 2, use_cuda = 0; end

if use_cuda
    X = transpose(dct(transpose(dct(X))));
else
    X = dct2(X);
end

end