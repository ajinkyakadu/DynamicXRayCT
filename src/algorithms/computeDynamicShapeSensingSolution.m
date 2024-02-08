function [V_iter, en_iter, t_iter] = computeDynamicShapeSensingSolution(A, y, n, numIter, use_cuda, ...
    options, init_sol, plot_iterates)
% GENERATEDYNAMICSHAPESENSINGSOLUTION generates a Dynamic Shape Sensing
% solution for the input data.
%
% Inputs:
%
% A : System matrix for tomography
% y : Observation data for tomography
% n : Dimension of the image
% numIter : Number of iterations for the optimization
% use_cuda : Flag for whether to use GPU or not
% options : Structure with optional parameters
% init_sol : Initial solution for the inversion
%
% Outputs:
%
% V_iter : Solution image
% phi_iter : Level-set function
% xpIter : Coefficient functions for the level-set
%
%
% See also minConf_SPG, levelSetDCT

% Define the default values for the options structure
if nargin < 6, options = []; end
if nargin < 8, plot_iterates = false; end

% Print the parameters and stats before the iterations start
fprintf('\n'); 
fprintf('=============================================================\n');
fprintf('<strong>                Dynamic Shape Sensing </strong>\n');
fprintf('Reconstruction with discrete and compressed sensing prior    \n');
fprintf('=============================================================\n');
disp('Parameters and Stats:');
disp(['Measurement Matrix        : ', num2str(size(A))]);
disp(['Measurement size          : ', num2str(size(y))]);
disp(['Dimension of the image    : ', num2str(n)]);
disp(['Number of iterations      : ', num2str(numIter)]);
disp(['Use CUDA                  : ', num2str(use_cuda)]);
fprintf('-------------------------------------------------------------\n');

% Display options fields
disp('DSS Options:');
fields = fieldnames(options);
for i = 1:numel(fields)
    fprintf('%25s : %s \n', fields{i}, num2str(options.(fields{i})));
end

if exist('init_sol', 'var')
    disp(['Initial solution size     : ', num2str(size(init_sol))]);
else
    disp('Initial solution not provided.');
end
if nargin < 8
    disp('Plot iterates: false');
else
    disp(['Plot iterates             : ', num2str(plot_iterates)]);
end

% Get the options for the optimization
kT       = getoptions(options, 'nBasis',    [0.1 0.1 0.1]);
kappa0   = getoptions(options, 'kappa',     0.1);
estimKp  = getoptions(options, 'estimKp',   1);
sigma    = getoptions(options, 'sigma',     0);
xEstIter = getoptions(options, 'xEstIter',  10);
kpEstIter= getoptions(options, 'kpEstIter', 10);
tau      = getoptions(options, 'tau',       Inf);

% generate the initial solution
if (~exist('init_sol', 'var') || isempty(init_sol))
    init_sol = generateInitialSolution(n, use_cuda);
end

% Initialize the coefficients
[x0, x_id] = initializeCoeff(init_sol - multithresh(gather(init_sol)), ...
    kT, use_cuda);

% Transfer the observation data to GPU if necessary
if use_cuda
    y = gpuArray(y);
end

% Define the tomography function handle
nY = 1; % norm(y(:))^2;
F  = @(x) misfitDynamicIP(x, A, y, nY, use_cuda);

% Define the projection function handle
funProj = @(x) projL1Ball(x, tau);

% Optimization options
optimOpt = struct();
optimOpt.maxIter    = xEstIter;
optimOpt.memory     = xEstIter;
optimOpt.progTol    = 1e-16;
optimOpt.testOpt    = 0;
optimOpt.optTol     = 1e-16;
optimOpt.print      = 0;
optimOpt.verbose    = 0;
optimOpt.curvilinear= 1;
optimOpt.bbType     = 1;

% Initializations
xIter = gather(x0);
kappa = kappa0;

% cast to GPU
% if use_cuda
%     xIter = gpuArray(xIter);
% end

% history structure
hist = struct();
hist.options = options;
hist.optimOpt= optimOpt;

DSSOpt = struct();

clock_cmp = tic;

% run multiple times with reducing boundary width
fprintf('-------------------------------------------------------------\n');
fprintf('computing reconstruction... \n');

for i=1:numIter
    
    % Define function handle
    DSSOpt.m0    = 0;
    DSSOpt.m1    = 1;
    DSSOpt.n     = n;
    DSSOpt.x_id  = gather(x_id);
    DSSOpt.sigma = sigma;
    DSSOpt.kappa = kappa;
    fh = @(x) levelSetDCT(x, F, use_cuda, DSSOpt);

    fprintf('iter:%2.0d/%2.0d kappa:%.4f tau=%8.2f/%8.2f, uval=%.2f | ', i, numIter, kappa, ...
        norm(xIter,1), tau, DSSOpt.m1);
    
    % minimize
    [xIter, ~] = minConf_SPG(fh, xIter, funProj, optimOpt);

    % Process intermediate images and update the optimization parameters
    [V_iter,phi_iter,~] = generateImage(xIter, n, x_id, kappa, use_cuda);

    % update the m1 value
    DSSOpt.m1 = updateIntensity(V_iter, A, y);
    
    if (i < numIter) && estimKp
        kappa = findKappa(kappa, phi_iter, A, y, kpEstIter, use_cuda);
    else
        kappa = 0.8*kappa;
    end
    
    % save history
    hist.f(i) = computeMainObjective(V_iter, A, y, use_cuda);
    % hist.g(i) = computeBoundaryPenalty(xIter, n, x_id, kappa, use_cuda);
    fprintf('f=%.4f \n', hist.f(i));
    
    % heaviside width tolerance condition
    if kappa < 1e-6, break; end
    
    % plot Figure
    if(plot_iterates)
        plotIntermediateSolution(V_iter)
    end
    
end

% Generate the final image
[V_iter, phi_iter, xpIter] = generateImage(xIter, n, x_id, 1e-12, use_cuda);

% Copy the data from the GPU to the MATLAB workspace
if use_cuda
    V_iter   = gather(V_iter);
    phi_iter = gather(phi_iter);
    xpIter   = gather(xpIter);
end

t_iter = toc(clock_cmp);
en_iter= hist.f;

fprintf('-------------------------------------------------------------\n');
fprintf('reconstruction size      : %d x %d x %d \n', n);
fprintf('Time taken               : %.4f \n', t_iter);
fprintf('=============================================================\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       Supporting Functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function uOpt = updateIntensity(V_iter, A, y)

Av = A*V_iter(:);
y  = y(:);

uOpt = (Av'*y)/norm(Av)^2;

end


function V = generateInitialSolution(n, use_cuda)
% generateInitialSolution generates an initial solution for a given n.
%
% Input:
% n - a vector of size 3 indicating the dimensions of the initial solution
% use_cuda - a boolean indicating whether to use GPU acceleration or not
%
% Output:
% V - the initial solution

% Initialize V with zeros
V = zeros(n, 'single');

% Generate a 3D grid of x, y, and z coordinates
[xx, yy, ~] = ndgrid(linspace(-1, 1, n(1)), linspace(-1, 1, n(2)), ...
    linspace(-1, 1, n(3)));

% Set all values of V within a circle of radius 1 to 1
V(xx.^2 + yy.^2 <= 1^2) = 1;

% If GPU acceleration is requested, convert V to a GPU array
if use_cuda
    V = gpuArray(V);
end

end



function [f] = computeMainObjective(V, A, y, use_cuda)
% compute the least-squares misfit
% Inputs:
% V - input spatiotemporal image
% A - tomography matrix
% y - tomography data
% use_cuda - flag for using GPU
%
% Outputs:
% f - objective value

% check if GPU acceleration is enabled
if use_cuda
    % compute the objective function with GPU acceleration
    f = gather(0.5 * norm(gpuArray(A * gather(V(:))) - y)^2);
else
    % compute the objective function without GPU acceleration
    f = 0.5 * norm(A * V(:) - y)^2;
end

end


function [penalty] = computeBoundaryPenalty(x, n, x_id, kappa, use_cuda)
%COMPUTEBOUNDARYPENALTY computes the boundary penalty
%
% Input:
%   x: current estimate of the image
%   n: size of the image
%   x_id: indices for the interior of the image
%   kappa: regularization parameter
%
% Output:
%   penalty: computed boundary penalty

% Initialize xv with the appropriate size
xv = zeros(prod(n), 1, 'single');

% Assign x to the interior of the image
if use_cuda
    xv(x_id) = gather(x);
end

% Transform xv using the DCT
z = dctThree(reshape(xv, n), use_cuda);

% Convert z to a vector
z = vec(z);

% Compute the Heaviside function and the Dirac delta function
options.epsi = kappa * (max(z) - min(z));
[~, delta] = heavi(z, options);

% Compute the boundary penalty
penalty = 0.5*norm(delta)^2;

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
phi = dctThree(xp, use_cuda);
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
inverseDCT = idctThree(V, use_cuda);

% Initialize the coefficients location matrix
coefficientLocation = zeros(size(inverseDCT));

% Calculate the number of coefficients to be considered
if length(kT) == length(imageSize)
    numCoeffs = floor(kT.*imageSize);
    coefficientLocation(1:numCoeffs(1),1:numCoeffs(2),1:numCoeffs(3)) = 1;
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


function [] = plotIntermediateSolution(V)
%PLOTINTERMEDIATESOLUTION Plot the intermediate solution.
% This function plots the intermediate solution for a video of size V.

% Get the number of frames in the video.
numFrames = size(V, 3);

% Choose 16 frames evenly spaced from the video to display.
selectedFrames = ceil(linspace(1, numFrames, 15));

% Open a new figure and display the selected frames as a mosaic.
figure(99);
out = imtile(V, 'Frames', selectedFrames, 'GridSize', [3 5]);
imshow(out);
pause(0.001);

end