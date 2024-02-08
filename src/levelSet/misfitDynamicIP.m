function [f, g] = misfitDynamicIP(X, A, Y, nY, use_cuda)
% MISFITDYNAMICIP - Computes the misfit function and gradient for dynamic 
% inverse problem
% 
% f(X) = \sum_{i} 0.5 * \| A_i X_i - Y_i \|_2^2 / nY
% g(X) = \sum_{i} A_i^T ( A_i X_i - Y_i) / nY
%
% [f,g] = MisfitDynamicIP(X, A, Y, nY, UseCuda) computes the misfit 
% function and gradient for the dynamic inverse problem, given the input 
% data X, measurement operator A, measured data Y, number of measurements 
% nY and flag to use CUDA or not UseCuda.
%
% Inputs:
%   X - The input data (n-dimensional array)
%   A - The measurement operator (matrix that maps from n-dimensional space
%       to m-dimensional space)
%   Y - The measured data (m-dimensional vector)
%   nY - The number of measurements (scalar)
%   use_cuda - A flag indicating whether to use CUDA (scalar, either 0 or 1)
%
% Outputs:
%   f - The misfit function value (scalar)
%   g - The gradient of the misfit function (n-dimensional array)

% Compute forward projections
if use_cuda
    Ax = gpuArray(A * gather(X));
else
    Ax = single(A * X);
end

% Calculate the residual
res = Ax - Y;

% Compute the misfit function value
f = 0.5 * norm(res)^2 / nY;

% Compute the gradient
if use_cuda
    Atr = gpuArray(A' * gather(res));
else
    Atr = single(A' * res);
end
g = Atr / nY;

end