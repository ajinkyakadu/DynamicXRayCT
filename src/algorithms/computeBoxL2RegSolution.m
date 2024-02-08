function [V_Temp, energy, t_Temp] = computeBoxL2RegSolution(A, y, n, beta, ...
    numIter, use_cuda, xMin, xMax)
% Generate Temporal Regularized Solution Function - This function 
% generates a temporal regularized solution by minimizing a misfit 
% function that includes a L2-regularization term on the difference of
% consequent time frames:
%
% V_Temp = argmin_x { \sum_t^T  \| A_t x_t - y_t \|_2^2 + \beta \| x_{t+1} - x_{t} \|_2^2 } 
%   subject to box constraints x_t \in [0,1]
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
if nargin < 4, beta     = 1e-3; end
if nargin < 5, numIter  = 100;  end
if nargin < 6, use_cuda = 0;    end 
if nargin < 7, xMin     = 0;    end
if nargin < 8, xMax     = 1;    end

fprintf('\n'); 
fprintf('=============================================================\n');
fprintf('<strong>                  BoxL2 Reconstruction </strong>\n');
fprintf('   Temporally regularized Least Squares with Box Constraints \n');
fprintf('=============================================================\n');
fprintf('MeasurementMatrix        : %d x %d \n', size(A));
fprintf('Measurement size         : %d x %d \n', size(y));
fprintf('Number of Iterations     : %d \n', numIter);
fprintf('Use GPU                  : %d \n', use_cuda);
fprintf('Motion reg parameter     : %d \n', beta);
fprintf('Bounding Box             : [%d, %d] \n', xMin, xMax);

% put the data on GPU if available
y = double(y);
if use_cuda, y = gpuArray(y); end

iterT    = 1;
numIterT = numIter;

% Define the function handle for L2 misfit and the bound projection function
normY= 1; % 0.5*norm(y(:))^2;
fh   = @(x) misfit(x, A, y, n, beta, normY, use_cuda);
fp   = @(x) boundProject(x, xMin, xMax);

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

fprintf('-------------------------------------------------------------\n');
fprintf('Optimization Algorithm   : %s \n', 'SPG');
fprintf('Memory                   : %d \n', optOptions.memory);
fprintf('Progress Tolerance       : %e \n', optOptions.progTol);
fprintf('Test OptimizationTol     : %d \n', optOptions.testOpt);
fprintf('Optimization Tolerance   : %d \n', optOptions.optTol);
fprintf('curvilinear              : %d \n', optOptions.curvilinear);
fprintf('Barzilai-Borwein Type    : %d \n', optOptions.bbType);
fprintf('-------------------------------------------------------------\n');
fprintf('computing reconstruction ... \n');

% Start progress indicator
progress('_start');
clock_cmp = tic;

% initial estimate
x0 = rescale(A'*gather(y));
if use_cuda, x0 = gpuArray(x0); end

% optimize using spectral gradient descent method
xopt = minConf_SPG(fh, x0, fp, optOptions);
energy = fh(xopt);

% gather the array from GPU
if use_cuda, xopt = gather(xopt); y = gather(y); end

% Reshape xr into the desired image size
V_Temp  = reshape(xopt,n);

% End progress display
progress('_end');
t_Temp = toc(clock_cmp);

fprintf('-------------------------------------------------------------\n');
fprintf('reconstruction size      : %d x %d x %d \n', n);
fprintf('Time taken               : %d \n', t_Temp);
fprintf('=============================================================\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supporting Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, g] = misfit(x, A, y, n, beta, normY, use_cuda)
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
K = single([1, -1]);
KT= single([-1, 1]);
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
f = 0.5 * (norm(res)^2 + beta * norm(Dx(:))^2)/normY;

% Compute the gradient
if use_cuda
    Atr = gpuArray(A' * gather(res));
else
    Atr = A'*res;
end

Dtr = conv2(Dx, KT, "full");
g   = Atr + beta * Dtr(:);
g   = g/normY;

f = gather(f);
g = gather(g);

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
