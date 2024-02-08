function [kappa] = findKappa(kappa0, phi, A, y, numIter, use_cuda)
% FINDKAPPA finds the optimal kappa (heaviside width) value for the 
% objective function (which is based on least-squares of tomographic
% measurements).
% 
% Inputs:
%   kappa0: Initial value for kappa.
%   phi: Level-set
%   A: Tomography matrix.
%   y: Tomography data.
%   numIter: number of iterations
%   use_cuda: flag for using GPU
% 
% Outputs:
%   kappa: Optimal kappa value.

% Define an anonymous function to compute the misfit function
fh = @(kappa) misfitKappa(kappa, phi, A, y, use_cuda);

% Set optimization options
options = optimset('MaxFunEvals', numIter, 'MaxIter', numIter, ...
    'Display', 'off');

% Minimize the misfit function using the given options
kappa = fminbnd(fh, 0, kappa0, options);

end

function [f] = misfitKappa(kappa, phi, A, y, use_cuda)
% MISFITKAPPA computes the misfit function for the given kappa value.
%   kappa: Current value of kappa.
%   phi: level-set function
%   A: Tomography matrix.
%   y: Tomography data.
%   use_cuda: flag to use GPU.
% 
% Output:
%   f: Value of the misfit function for the given kappa.

% Vectorize the DCT coefficients
phi = phi(:);

% Compute the Heaviside function
hopt.epsi = kappa * (max(phi) - min(phi));
[h] = heavi(phi, hopt);

% compute the residual
if use_cuda
    res = gpuArray(A*gather(h)) - y;
else
    res = A*h - y;
end

% Compute the misfit function & return the value of the misfit function
if use_cuda
    f = double(gather(norm(res)^2));
else
    f = double(norm(res)^2);
end

end


