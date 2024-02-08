function [V_Temp] = generateBoxL2RegSolution(A, y, n, lambda, ...
    numIter, use_cuda)
% Generate Temporal Regularized Solution Function - This function 
% generates a temporal regularized solution by minimizing a misfit 
% function that includes a L2-regularization term on the difference of
% consequent time frames.
%
% Inputs:
% A - A matrix
% y - A vector
% n - The size of the desired output image
% lambda (optional) - The regularization parameter (default: 1e-3)
% numIter (optional) - The maximum number of iterations for the 
%   optimization (default: 100)
% use_cuda (optional) - binary value that indicates whether to use GPU or 
%   not (default value is 0)
%
% Outputs:
% V_Temp - The generated temporal regularized solution
global iterT numIterT

% Check if the lambda and numIter inputs were provided
if nargin < 4, lambda   = 1e-3; end
if nargin < 5, numIter  = 100; end
if nargin < 6, use_cuda = 0; end 

fprintf('Generating Box-L2 regularized solution...\n');

% put the data on GPU if available
y = double(y);
if use_cuda, y = gpuArray(y); end

% Start progress indicator
progress('_start');
iterT    = 1;
numIterT = numIter;

% Define the function handle for L2 misfit and the bound projection function
normY= 1; % 0.5*norm(y(:))^2;
fh   = @(x) misfit(x, A, y, n, lambda, normY, use_cuda);
fp   = @(x) boundProject(x, 0, 1);

% optimization options
optOptions = [];
optOptions.maxIter      = numIter;
optOptions.memory       = 10;
optOptions.progTol      = 1e-16;
optOptions.testOpt      = 1;
optOptions.optTol       = 1e-16;
optOptions.print        = 0;
optOptions.verbose      = 0;
optOptions.curvilinear  = 1;
optOptions.bbType       = 1;


% initial estimate
x0 = rescale(A'*gather(y));
if use_cuda, x0 = gpuArray(x0); end

% optimize using spectral gradient descent method
xopt = minConf_SPG(fh, x0, fp, optOptions);

% gather the array from GPU
if use_cuda, xopt = gather(xopt); y = gather(y); end

% Reshape xr into the desired image size
V_Temp  = reshape(xopt,n);

% End progress display
progress('_end');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supporting Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, g] = misfit(x, A, y, n, lambda, normY, use_cuda)
% Misfit function - Computes the objective function and gradient
%
% Inputs:
% x      - Input vector
% A      - Matrix A
% y      - Target vector
% n      - Size of the input
% lambda - Regularization parameter
% use_cuda - binary value that indicates whether to use GPU or not 
%
% Outputs:
% f - Objective function value
% g - Gradient of the objective function
global iterT numIterT

% define forward and backward kernel for temporal regularization
K = single([1, -1] / sqrt(2));
KT= single([-1, 1] / sqrt(2));
if use_cuda, K = gpuArray(K); KT = gpuArray(KT); end

% Calculate the difference along the temporal direction
Dx = conv2(reshape(x, [], n(3)), K, "valid");

% Calculate the product of matrix A and x
if use_cuda
    Ax = gpuArray(A * gather(x));
else
    Ax = A*x;
end

% Calculate the residual
res = Ax - y;

% Compute the objective function value
f = 0.5 * (norm(res)^2 + lambda * norm(Dx(:))^2)/normY;

% Compute the gradient
if use_cuda
    Atr = gpuArray(A' * gather(res));
else
    Atr = A'*res;
end

Dtr = conv2(Dx, KT, "full");
g   = Atr + lambda * Dtr(:);
g   = g/normY;

iterT = iterT + 1;
progress(iterT, numIterT)

end

function [x] = boundProject(x, LB, UB)
% Bounds the elements of x to the lower bound LB and upper bound UB
%
% Inputs:
% x  - Input vector
% LB - Lower bound
% UB - Upper bound
%
% Outputs:
% x - Bounded vector

% Clip the values to the lower and upper bounds
x(x < LB) = LB;
x(x > UB) = UB;

end
