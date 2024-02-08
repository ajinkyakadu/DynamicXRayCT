function [V_Temp,info] = computeBoxTVL2RegSolution(A, y, n, alpha, beta, numIter)
% Generate Temporal Regularized Solution Function - This function 
% generates a temporal regularized solution by solving this optimization
% problem:
%
% V_Temp = argmin_x { \sum_t^T  \| A_t x_t - y_t \|_2^2 + \alpha \|\nabla x_t \|_1 + \beta \|x_{t+1} - x_{t} \|_2^2 } 
%   subject to box constraints x_t \in [0,1]
%
% Inputs:
%   A - A matrix
%   y - A vector
%   n - The size of the desired output image
%   alpha (optional) - The regularization parameter (default: 1e-3)
%   beta (optional)  - The regularization parameter (default: 1e-3)
%   numIter (optional) - The maximum number of iterations for the 
%   optimization (default: 100)
%   use_cuda (optional) - binary value that indicates whether to use GPU or 
%   not (default value is 0)
%
% Outputs:
% V_Temp - The generated temporal regularized solution

% Check if the lambda and numIter inputs were provided
if nargin < 4, alpha   = 1*10^-3; end
if nargin < 5, beta    = 1*10^-3; end
if nargin < 6, numIter = 100; end 

fprintf('Generating Box-TV-L2 regularized solution...\n');

% put the data on GPU if available
y = single(y);

for t = 1:n(3)
    Af{t}   = @(x) A{t} * x(:);
    ATf{t}  = @(x) reshape((A{t}' * x), [n(1),n(2)]);
end

dynA  = @(x) dynamicOperator(Af, x);
dynAT = @(f) dynamicOperator(ATf, f);



optimOpt               = [];
optimOpt.OFtype        = 'linear';
optimOpt.TVtypeU       = 'isotropic';
optimOpt.constraint    = 'range';
optimOpt.conRange      = [0, 1];
optimOpt.computeEnergy = true;
optimOpt.szU           = n;
optimOpt.useATF        = false;
optimOpt.p             = 2;
optimOpt.rho           = 1;
optimOpt.maxEval       = numIter;
optimOpt.rhoAdaptation = true;
optimOpt.rhoAdaptK     = 25;
optimOpt.overRelaxPara = 1.8;
optimOpt.dataCast      = 'double';
optimOpt.output        = true;
optimOpt.regCon = 1;
optimOpt.lsSolverPara  = [];
optimOpt.lsSolverPara.lsSolver      = 'CG';
optimOpt.lsSolverPara.stopCriterion = 'progRelRes';
optimOpt.lsSolverPara.progMode      = 'poly';
optimOpt.lsSolverPara.tol           = 10^-3;
optimOpt.lsSolverPara.tolDecExp     = 1.5;
optimOpt.lsSolverPara.minIter       = 10;
optimOpt.lsSolverPara.maxIter       = numIter/2;

v              = zeros([n(1), n(2), 2, n(3)-1], 'like', y);
% due to different conventions, we need to divide alpha by 2 and sqrt
% beta
[V_Temp, info] = TVOF_Deblurring_ADMM(dynA, dynAT, y, v, alpha/2, sqrt(beta), optimOpt);

end
