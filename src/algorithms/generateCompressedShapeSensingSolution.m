function [V_frame] = generateCompressedShapeSensingSolution(A, y, n, nB, numIter, use_cuda, options)
%generateFbyFSolution generates a frame-by-frame solution to following 
%                     optimization problem.
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
if nargin < 4, nB       = 1e-3; end
if nargin < 5, numIter  = 100;  end
if nargin < 6, use_cuda = 0;    end
if nargin < 7, options  = [];   end

% Display progress message
fprintf('Generating Compressed Shape Sensing solution (nB = %d) ...\n', nB);

% Define number of frames
batches = createBatches(n(3),nB);
nFrames = length(batches);
V_frame = zeros(n,'single');

% Start progress indicator
progress('_start')

% Loop through each frame
for i=1:nFrames
    
    % Extract the current frame
    currentY = y(:,batches{i});
    currentA = A(batches{i});
    
    % Concatenate the tomography matrices
    AI = currentA{1};
    for j=2:length(currentA)
        AI = [AI;currentA{j}];
    end
    
    % Vectorize the current frame & move to GPU (if required)
    currentY = single(vec(currentY));
    if use_cuda
        currentY = gpuArray(currentY);
    end

    % optimize here
    xResult = solveStaticCSSProblem(AI, currentY, [n(1) n(2)], numIter, use_cuda, options);
    
    % Reshape x into the desired image size
    V_frame(:,:,batches{i}) = repmat(reshape(xResult, [n(1) n(2)]), [1 1 nB]);

    % Display progress
    progress(i, nFrames);
end

% End progress display
progress('_end');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Supporting Functions
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


function [xResult] = solveStaticCSSProblem(AI, currentY, n, numIter, use_cuda, options)

% Get the options for the level-set
kT       = getoptions(options, 'nBasis',    [0.1 0.1]);
kappa0   = getoptions(options, 'kappa',     0.1);
xEstIter = getoptions(options, 'xEstIter',  20);
tau      = getoptions(options, 'tau',       Inf);
uval     = getoptions(options, 'uval',      1);
% Initialize x0
I0 = generateInitialSolution(n, use_cuda);
[x0, x_id] = initializeCoeff(I0, kT, use_cuda);

% Define the function handle for A transpose times A
nY = 0.5*norm(currentY(:))^2;
misfitFunc  = @(x,kappa) misfitCSS(x, x_id, AI, currentY, n, nY, kappa, uval, use_cuda);
projectFunc = @(x) projL1Ball(x, tau);

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

kappaI = kappa0;
xIter  = x0;

for i=1:floor(numIter/xEstIter)

    fh = @(x) misfitFunc(x, kappaI);

    % minimize
    [xIter, ~] = minConf_SPG(fh, xIter, projectFunc, optimOpt);

    % update kappa
    kappaI = 0.9*kappaI;
    
    % heaviside width tolerance condition
    if kappaI < 1e-3, break; end

end

% Generate the final image
xResult = generateImage(xIter, n, x_id, 1e-16, use_cuda);

if use_cuda, xResult = gather(xResult); end

end

function [f, g] = misfitCSS(x, x_id, A, y, n, nY, kappa, m1, use_cuda)
% DESCRIPTION:
%   This function implements the computation of the objective function and 
%   its gradient for the Level Set DCT problem.
%
% INPUTS:
%   x       : Vector of DCT coefficients
%   funObj  : Function handle to compute the objective and gradient of the
%             underlying image
%   use_cuda: Boolean flag to enable GPU acceleration
%   options : Structure containing options to configure the computation
%
% OUTPUTS:
%   f   : Value of the objective function
%   g   : Gradient of the objective function wrt DCT coefficients
%   m   : Computed intensity scaled image
%   d   : Dirac-delta function computed from the DCT coefficients


% Convert to the image domain
xv = zeros(prod(n), 1, 'single');
if use_cuda, xv = gpuArray(xv); end
xv(x_id) = x;

z = dctTwo(reshape(xv, n), use_cuda);
z = vec(z);

% Heaviside width
hopt.epsi = kappa * (max(z) - min(z));

% Function and gradient
[h, d] = heavi(z, hopt);


% Scale the image by intensity
m = h * m1;

% Compute the objective and gradient wrt image
if use_cuda
    res = gpuArray(A*gather(m)) - y;
else
    res = A*m - y;
end

f   = 0.5*norm(res)^2/nY;

if use_cuda
    g0 = gpuArray(A'*gather(res))/nY;
else
    g0  = A'*res/nY;
end
g0  = m1 * (g0 .* d);

% Compute gradient wrt DCT coefficients
g0 = vec(idctTwo(reshape(g0, n), use_cuda));
g  = g0(x_id);

end


function V = generateInitialSolution(n, use_cuda)
% generateInitialSolution generates an initial solution for a given n.
%
% Input:
% n - a vector of size 2 indicating the dimensions of the initial solution
% use_cuda - a boolean indicating whether to use GPU acceleration or not
%
% Output:
% V - the initial solution

% Initialize V with zeros
V = zeros(n, 'single');

% Generate a 3D grid of x, y, and z coordinates
[xx, yy] = ndgrid(linspace(-1, 1, n(1)), linspace(-1, 1, n(2)));

% Set all values of V within a circle of radius 1 to 1
V(xx.^2 + yy.^2 <= 1^2) = 1;

% If GPU acceleration is requested, convert V to a GPU array
if use_cuda
    V = gpuArray(V);
end

end

function [V, phi, xp] = generateImage(x, n, x_id, kappa, use_cuda)
% Generates an image, a phase map and an intermediate result
%
% Inputs:
%   x - an array representing the input data
%   n - the dimensions of the resulting image
%   x_id - the index of the input data in the full array
%   kappa - a scalar value used in the calculation
%
% Outputs:
%   V - the resulting image
%   phi - the phase map
%   xp - the intermediate result

% Initialize xp to a zero array with the specified dimensions
xp = zeros(prod(n), 1, 'single');
if use_cuda
    xp = gpuArray(xp);
end

% Set the values of x at the specified indices to the input data
xp(x_id) = x;

% Reshape xp to the specified dimensions
xp = reshape(xp, n);

% Compute the DCT
phi = dctTwo(xp, use_cuda);
phi = vec(phi);

% Apply the Heaviside function
hopt.epsi = kappa * (max(phi) - min(phi));
V = heavi(phi, hopt);

% Reshape the resulting image to the specified dimensions
V = reshape(V, n);

end


function [x,x_id] = initializeCoeff(V, kT, use_cuda)
% INPUTS:
% V - A 3D array representing the input image
% kT - The fraction of the coefficients to be considered
%
% OUTPUTS:
% x0 - A vector containing the restricted IDCT coefficients
% x_id - A vector containing the indices of the restricted coefficients

% Get the size of the image
imageSize = size(V);

% Compute the inverse DCT of the image
inverseDCT = idctTwo(V, use_cuda);

% Initialize the coefficients location matrix
coefficientLocation = zeros(size(inverseDCT));

% Calculate the number of coefficients to be considered
if length(kT) == length(imageSize)
    numCoeffs = floor(kT.*imageSize);
    coefficientLocation(1:numCoeffs(1),1:numCoeffs(2)) = 1;
else
    numCoeffs = floor(kTprod(imageSize));
    [~, xda_sort_idx] = sort(abs(inverseDCT(:)), 'descend');
    coefficientLocation(xda_sort_idx(1:numCoeffs)) = 1;
end

% Restrict the inverse DCT coefficients
xr = inverseDCT.*coefficientLocation;

% Vectorize the coefficients
xr = xr(:);

% Find the indices and values of the restricted coefficients
x_id = find(xr);
x    = xr(x_id);

end
