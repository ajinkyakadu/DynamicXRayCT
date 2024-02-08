function [V_frame, misfit_frame, t_frame] = computeStaticSolution(A, y, n, nB, numIter, use_cuda, use_FL, alpha_TV)
%computeStaticSolution generates a frame-by-frame solution to following
%                     optimization problem.
%
%   computeStaticSolution computes a bin-by-bin reconstruction solving 
%   \| A x - f \|_2^2 + alpha TV(x)
%   and allows for the inclusion of box constraints
%
% Summary: The function solves the optimization problem by processing nB
% frames at a time, where nB is the number of frames specified by the user
% or the default value of 1e-3. The optimization problem is solved using
% the preconditioned conjugate gradient method with the minConf_SPG
% function. The result is stored as a 3D matrix with the size specified
% by the n input.
%
% INPUTS:
%   - A: a cell array that contains a matrix
%   - y: a matrix that contains the measurement data
%   - n: a 3-element vector that defines the size of the result
%   - nB: scalar that defines the number of frames to process at a time
%           (default value is 1e-3)
%   - numIter: integer that defines the maximum number of iterations for
%           optimization (default value is 100)
%   - use_cuda: binary value that indicates whether to use GPU or not
%           (default value is 0)
%
% OUTPUTS:
%   - V_frame: a 3D matrix that contains the result of the optimization problem

% Check for optional input arguments
if nargin < 4, nB       = 1;    end
if nargin < 5, numIter  = 100;  end
if nargin < 6, use_cuda = 0;    end
if nargin < 7, use_FL   = 0;    end
if nargin < 8, alpha_TV = 0;    end


if (alpha_TV > 0)
    use_FL = true;
end

% Define number of frames
batches       = createBatches(n(3),nB);
nFrames       = length(batches);
V_frame       = zeros(n,'single');
misfit_frame  = zeros(1,nFrames);

% Display progress message
fprintf('\n'); 
fprintf('=============================================================\n');
fprintf('<strong>                Static Binned Reconstruction</strong>\n');
fprintf('                     (Bin size: %d)\n', nB);
fprintf('=============================================================\n');
fprintf('MeasurementMatrix/Frame  : %d x %d \n', size(A{1}));
fprintf('Number of total frames   : %d \n', length(A));
fprintf('Measurement size         : %d x %d \n', size(y));
fprintf('Batch size               : %d \n', nB);
fprintf('Number of Iterations     : %d \n', numIter);
fprintf('Use GPU                  : %d \n', use_cuda);
fprintf('TV reg parameter         : %d \n', alpha_TV);
fprintf('number of Bins           : %d \n', nFrames);

fprintf('-------------------------------------------------------------\n');
fprintf('computing reconstruction for each bin... \n');
% Start progress indicator
progress('_start')


% Loop through each frame
clock_cmp = tic;
for i=1:nFrames
    
    % Extract the current frame
    currentY = y(:,batches{i});
    currentA = A(batches{i});
    
    % Concatenate the tomography matrices
    AI = currentA{1};
    for j=2:length(currentA)
        AI = [AI;currentA{j}];
    end
    
    % Vectorize the current frame & Move data to GPU (if required)
    currentY = vec(currentY);
    currentY = single(currentY);
    if use_cuda
        currentY = gpuArray(currentY);
    end
    
    % optimize here
    [xResult, misfit_frame(i)] = solveStaticProblem(AI, currentY, [n(1) n(2)], ...
        numIter, use_cuda, use_FL, alpha_TV);
    
    % Reshape x into the desired image size
    V_frame(:,:,batches{i}) = repmat(reshape(gather(xResult), [n(1) n(2)]), [1 1 nB]);
    
    % Display progress
    progress(i, nFrames);
    
end

% End progress display
progress('_end');
t_frame = toc(clock_cmp);

fprintf('-------------------------------------------------------------\n');
fprintf('reconstruction size      : %d x %d x %d \n', n);
fprintf('Time taken               : %.4f \n', t_frame);
fprintf('=============================================================\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Supporting Functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [batches] = createBatches(N, M)

% Compute the number of batches and the size of the last batch
num_batches = ceil(N/M);
last_batch_size = mod(N,M);

% Create a cell array to store the batches
batches = cell(num_batches,1);

% Loop over the batches and fill them with the elements of the linear array
for i = 1:num_batches
    if i == num_batches && last_batch_size ~= 0 % Last batch with overlap
        batches{i} = (N-M+1):N;
    else % Other batches without overlap
        batches{i} = (i-1)*M+1 : min(i*M, N);
    end
end

% % Print the result
% for i = 1:num_batches
%     fprintf('Batch %d: %s\n', i, mat2str(batches{i}));
% end

end


function [xResult, misfitResult] = solveStaticProblem(AI, currentY, n, numIter, use_cuda, use_FL, alpha_TV)

% Define the function handle for A transpose times A
nY = 0.5*norm(currentY(:))^2;
misfitFunc  = @(x) misfit(x, AI, currentY, nY, use_cuda);
projectFunc = @(x) boundProject(x, 0, 1);

% Initialize x0
x0 = zeros(prod(n),1,'single');


% Solve for x using the preconditioned conjugate gradient method
if(use_FL)
    
    optimOpt = [];
    optimOpt.maxIter         = numIter;
    optimOpt.output          = false;
    optimOpt.displayWarnings = false;
    Proxy = @(x, nu, s, update) projBoxConstraints(x,  'range', [0,1]);
    
    if(alpha_TV == 0)
        
        optimOpt.J               = @(res) 0;
        optimOpt.Jx              = 0;
        optimOpt.fastGradient    = true;
        
        % compute lipschitz constant
        lipConstant = powerIteration(@(x) AI' * (AI * x), [prod(n), 1], 10^-10, 1, 0);
        nu                      = 1/lipConstant; % we normalized A
        
        % accelerated projected gradient algorithm
        xResult = ProxGradLinearLeastSquares(...
            @(x) AI*x,  @(x) AI'*x, gather(currentY), @(x) x, Proxy, nu, x0, ...
            zeros(size(currentY)), AI'*gather(currentY), optimOpt);
        
    else
        
        optimOpt.acceleration = false;
        optimOpt.constraint   = 'range';
        optimOpt.conRange     = [0,1];
        optimOpt.LipPowerIterTol = 10^-10;
        % we need to divide alpha by 2, due to different conventions
        [xResult, ~,~,info]  = TV_Deblurring(@(x) AI*x(:), @(x) reshape(AI'*x, n(1:2)), gather(currentY), alpha_TV/2, optimOpt);
        xResult  = xResult(:);
        
    end
    
else
    
    if use_cuda
        x0 = gpuArray(x0);
    end
    
    
    % Optimization options
    optimOpt = [];
    optimOpt.maxIter = numIter;
    optimOpt.memory  = 10;
    optimOpt.progTol = 1e-16;
    optimOpt.testOpt = 0;
    optimOpt.optTol  = 1e-16;
    optimOpt.print   = 0;
    optimOpt.verbose = 0;
    optimOpt.curvilinear = 1;
    optimOpt.bbType  = 1;
    
    xResult = minConf_SPG(misfitFunc, x0, projectFunc, optimOpt);
    
end

misfitResult = misfitFunc(xResult);

end


function [f, g] = misfit(x, A, y, nY, use_cuda)
% The function computes the residual and gradient of the misfit function
% between the predicted value and the observed data.
%
% Inputs:
% x - The current estimate of the parameters
% A - The design matrix
% y - The observed data
% nY - The number of observations
% use_cuda - Boolean flag indicating whether to use the GPU for computation

% Outputs:
% f - The value of the misfit function
% g - The gradient of the misfit function

% Check if CUDA should be used
if use_cuda
    Ax = gpuArray(A * gather(x));
else
    Ax = A*x;
end

% Compute the value of the misfit function
res = Ax - y;

% Calculate function value
f = gather(0.5 * norm(res)^2 / nY);

% Calculate gradient
if use_cuda
    Atr = gpuArray(A' * gather(res));
else
    Atr = A'*res;
end

% Compute the gradient of the misfit function
g   = gather(Atr / nY);

end


function [x] = boundProject(x, LB, UB)
% The function projects the input vector to the lower and upper bounds.
% Inputs:
% x - The input vector
% LB - The lower bound
% UB - The upper bound
%
% Outputs:
% x - The projected vector

% Clip values to lower and upper bounds
x(x < LB) = LB;
x(x > UB) = UB;

end