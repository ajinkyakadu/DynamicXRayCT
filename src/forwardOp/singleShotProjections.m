function [A, AC, y, thVec] = singleShotProjections(V_true, dTheta, noiseL, use_cuda)
% SINGLESHOTPROJECTIONS Generates the measurements from the true image 
%  using parallel-beam projections and adds noise to the measurements.
%
% INPUT: 
%   - V_true: The true temporal image with dimensions nx x ny x nT.
%   - dTheta: The measurement angle increment in degrees.
%   - noiseL: The function handle for adding noise to the measurements.
%
% OUTPUT:
%   - A:      The measurement matrix of size n_y x n_x*nT.
%   - AC:     The cell array of measurement matrices for each frame.
%   - y:      The noisy measurements of size n_z x nT.
%   - thVec:  The measurement angles in degrees.

if nargin < 4, use_cuda = 0; end


% Get the number of frames and the size of the images
nT = size(V_true,3);
n  = [size(V_true,1) size(V_true,2)];

% Display the start of generating measurements
fprintf('\n'); 
fprintf('=============================================================\n');
fprintf('<strong>         Generating Tomographic Measurements</strong>\n');
fprintf('                (Parallel-beam X-ray) \n');
fprintf('=============================================================\n');
fprintf('Image Size        : %d x %d \n', n);
fprintf('Timestamps        : %d \n', nT)

% Calculate the measurement angles
thVec = 0:dTheta:dTheta*(nT-1);
fprintf('Projections angles: %d : %d : %d (in degrees) \n', thVec(1), dTheta, thVec(end));

% Initialize an empty array to store the measurements
y = [];
AC= [];

% Set the volume geometry
vol_geom = astra_create_vol_geom(n(1), n(2));

% Generate projection data for each angle
for i = 1:nT
    theta = thVec(i);

    % Set the projection geometry for the current angle
    proj_geom = astra_create_proj_geom('parallel', 1, ceil(sqrt(2)*n(1)), theta*(pi/180));
    
    % set up Tomo SPOT operator
    if (use_cuda)
        kernelInv = 'cuda';
    else
        kernelInv = 'linear';
    end

    kernelFwd = 'strip';
    Ainv = opTomo(kernelInv, proj_geom, vol_geom);
    Afwd = opTomo(kernelFwd, proj_geom, vol_geom);
    
    % merge for inversion
    AC{i}     = Ainv;
    
    % Generate the measurement for the current angle and add noise
    y0 = Afwd*vec(V_true(:,:,i));
    ry = randn(size(y0));
    y(:,i) = y0 + noiseL*(norm(y0)/norm(ry))*ry;

    if i==1
        fprintf('Geometry          : %s \n', 'parallel');
        fprintf('detectors         : %d \n', ceil(sqrt(2)*n(1)));
        fprintf('Forward Kernel    : %s \n', kernelFwd);
        fprintf('Inverse Kernel    : %s \n', kernelInv);
        fprintf('Noise Type        : %s \n', 'Additive Gaussian');
        fprintf('Noise Strength    : %d \n', noiseL);
    end

end

% Form the block diagonal measurement matrix
A = blkdiag(AC{:});
y = single(y);

% Display the size of the measurement matrix and the number of measurements per frame
fprintf('Measurement matrix: %d x %d \n', size(A));
fprintf('measurements/frame: %d (single-shot)\n', size(y,1));
fprintf('=============================================================\n');

end

