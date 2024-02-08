function [A, Af, data, angles] = getRealTomographicMeasurements(stack, nT, angDif, use_cuda)
% getRealTomographicMeasurements - Obtains the real tomographic 
% measurements from the input data stack.
% Inputs:
%   stack   - A struct containing the following fields:
%   .rec    - 3D reconstruction volume
%   .angles - 2D array of projection angles
%   .data   - 3D array of projection data
%   .decPixelSz - Scalar value of detector pixel size
%   .detectorSz - Scalar value of detector size
%   .SOD    - Scalar value of source to origin distance
%   .SDD    - Scalar value of source to detector distance
%   nT      - Scalar value of the number of time steps
%   angDif  - Scalar value of the angle difference between two consecutive
%               projections
% 
% Outputs:
%   A     - Forward operator in a block diagonal form
%   Af    - Cell array of forward operators, one for each time step
%   data  - 2D array of projection data
%   angles- 1D array of projection angles
%
% Example:
% [A, Af, data, angles] = getRealTomographicMeasurements(stack, nT, angDif);

if nargin < 4, use_cuda = 0; end


% convert to single
stack.rec = im2single(stack.rec);
stack.data= im2single(stack.data);

% Obtain the size of the reconstruction
n = size(stack.rec);

% Obtain the projection angles
anglesJ = stack.angles;

% Define the volume geometry
volGeo = astra_create_vol_geom(n(1), n(2));

% Initialize the data and angle arrays
data   = [];
angles = zeros(1, nT);

% Loop over all time steps
for i = 1:nT
    % Compute the index of the projection angle for the current time step
    angleIdx = floor(rem(size(anglesJ, 1) * (angDif * i / 360), size(anglesJ, 1))) + 1;

    % Obtain the projection angle and data for the current time step
    angleI = anglesJ(angleIdx, i);
    data(:, i) = vec(stack.data(angleIdx, :, i));
    angles(i) = angleI;

    % Define the projection geometry
    projGeo = astra_create_proj_geom('fanflat', stack.decPixelSz, ...
        stack.detectorSz, angleI, stack.SOD, stack.SDD - stack.SOD);
    projGeo = astra_geom_2vec(projGeo);

    % Obtain the forward operator for the current time step
    if(use_cuda)
        Af{i} = opTomo('cuda', projGeo, volGeo);
    else
        Af{i} = opTomo('line_fanflat', projGeo, volGeo);
    end

end


% Obtain the forward operator in block diagonal form
A = blkdiag(Af{:});

% rescale data
data    = rescale(data);
fwdProj = A*vec(stack.rec(:,:,1:nT));
alpha   = (fwdProj'*data(:))/norm(data(:))^2;
data    = alpha*data;
data    = single(data);

end