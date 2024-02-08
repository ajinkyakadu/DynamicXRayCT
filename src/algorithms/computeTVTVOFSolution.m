function [V_Temp, U_Temp, info] = computeTVTVOFSolution(A, y, n, alpha, beta, gamma, ...
    numIterOuter, numIterImage, numIterMotion)
% Generates Motion regularized Tomographic solution - This function
% generates a motion regularized solution by minimizing a misfit
% function that includes a L2-regularization term on optical flow
% constraints (that captures the motion using velocity field) along with
% the total-variation norm on spatial images as well as total-variation
% norm on the velocity vector fields:
%
% (V_Temp, U_motion) = argmin_(x,u) { \sum_t^T  \| A_t x_t - y_t \|_2^2
%   + \alpha \|\nabla x_t \|_1 
%   + \beta \| x_{t+1} - x_t + (\nabla u_{t+1}) \cdot x_t \|_2^2 
%   + \gamma \sum_i^d \|\nabla u_{i,t} \|_1}
%   subject to box constraints x_t \in [0,1]
%
% Inputs:
%   A - A matrix
%   y - A vector
%   n - The size of the desired output image
%   alpha (optional) - The regularization parameter (default: 1e-3)
%   beta (optional) - The regularization parameter (default: 1e-3)
%   gamma (optional) - The regularization parameter (default: 1e-3)
%   numIter (optional) - The maximum number of iterations for the
%   optimization (default: 100)
%   use_cuda (optional) - binary value that indicates whether to use GPU or
%   not (default value is 0)
%
% Outputs:
% V_Temp - The generated temporal regularized solution

% Check if the lambda and numIter inputs were provided
if nargin < 4, alpha   = 2*10^-3; end
if nargin < 5, beta    = 2*10^-3; end
if nargin < 6, gamma    = 2*10^-3; end
if nargin < 7, numIterOuter = 4; end
if nargin < 8, numIterImage = 1000; end
if nargin < 9, numIterMotion = 1000; end

fprintf('Generating Box-TV-TV-OF regularized solution...\n');

% put the data on GPU if available
y = single(y);

for t = 1:n(3)
    Af{t}   = @(x) A{t} * x(:);
    ATf{t}  = @(x) reshape((A{t}' * x), [n(1),n(2)]);
end

dynA  = @(x) dynamicOperator(Af, x);
dynAT = @(f) dynamicOperator(ATf, f);


% set up regularization and optimization
optimOpt = [];
optimOpt.maxIter         = numIterOuter;
optimOpt.constraint      = 'range';
optimOpt.OFtype          = 'linear';
optimOpt.TVtypeU         = 'isotropic';
optimOpt.constraint      = 'range';
optimOpt.conRange        = [0, 1];
optimOpt.dimU            = n;
optimOpt.useATF          = false;
optimOpt.displayWarnings = true;

% general parameter of the alternation
optimOpt.hardCompRestr   = false;
optimOpt.output          = true;
optimOpt.monotoneEnergy  = false;
optimOpt.symmetricOF     = false;
optimOpt.p               = 2;

% x para (called u in Felix code)
optimOpt.uOpt            = [];
optimOpt.uOpt.solver     = 'ADMM';
optimOpt.uOpt.algo       = [];
optimOpt.uOpt.algo.rho           = 10^0;
optimOpt.uOpt.algo.maxEval       = numIterImage;
optimOpt.uOpt.algo.rhoAdaptation = true;
optimOpt.uOpt.algo.rhoAdaptK     = 25;
optimOpt.uOpt.algo.overRelaxPara = 1.8;
optimOpt.uOpt.algo.output        = false;
optimOpt.uOpt.algo.lsSolverPara  = [];
optimOpt.uOpt.algo.lsSolverPara.lsSolver      = 'CG';
optimOpt.uOpt.algo.lsSolverPara.stopCriterion = 'progRelRes';
optimOpt.uOpt.algo.lsSolverPara.progMode      = 'poly';
optimOpt.uOpt.algo.lsSolverPara.tol           = 10^-3;
optimOpt.uOpt.algo.lsSolverPara.tolDecExp     = 1.5;
optimOpt.uOpt.algo.lsSolverPara.minIter       = 3;
optimOpt.uOpt.algo.lsSolverPara.maxIter       = numIterImage;

% u para (called v in Felix code)
optimOpt.vOpt                  = [];
optimOpt.vOpt.solver           = 'ADMM';
optimOpt.vOpt.parallel         = true;
optimOpt.vOpt.nWorkerPool      = min(n(3)-1,maxNumCompThreads());
optimOpt.vOpt.algo             = [];
optimOpt.vOpt.algo.output           = false;
optimOpt.vOpt.algo.maxEval          = numIterMotion;
optimOpt.vOpt.algo.rho              = 10^0;
optimOpt.vOpt.algo.rhoAdaptation    = false;
optimOpt.vOpt.algo.rhoAdaptK        = 25;
optimOpt.vOpt.algo.overRelaxPara    = 1.8;
optimOpt.vOpt.algo.returnAlgoVar    = true;
optimOpt.vOpt.algo.lsSolverPara               = [];
optimOpt.vOpt.algo.lsSolverPara.lsSolver      = 'AMG-CG';
optimOpt.vOpt.algo.lsSolverPara.stopCriterion = 'progRelRes';
optimOpt.vOpt.algo.lsSolverPara.progMode      = 'poly';
optimOpt.vOpt.algo.lsSolverPara.tol           = 10^-4;
optimOpt.vOpt.algo.lsSolverPara.tolDecExp     = 1.5;
optimOpt.vOpt.algo.lsSolverPara.minIter       = 3;
optimOpt.vOpt.algo.lsSolverPara.maxIter       = 100;

% due to different conventions, we need to divide alpha and gamma by 2 and sqrt
% beta. Also note the different order
[V_Temp, U_Temp, info] = TVTVOF_Deblurring(dynA, dynAT, y, alpha/2, gamma/2, sqrt(beta), optimOpt);

V_Temp = single(V_Temp);
U_Temp = single(U_Temp);
end
