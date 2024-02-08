function X = idctTwo(X, use_cuda)
% IDCTThree - Perform Inverse Discrete Cosine Transform along two 
% dimensions
%
% Syntax: X = idctTwo(X)
%
% Inputs:
%   X - Input data to be transformed
%
% Outputs:
%   X - Transformed data
%
% The function performs the inverse discrete cosine transform along two 
% dimensions. The function also casts the input to the GPU to ensure 
% computation takes place on the GPU.

if nargin < 2, use_cuda = 0; end

if use_cuda
    X = idct(transpose(idct(transpose(X))));
else
    X = idct2(X);
end

end