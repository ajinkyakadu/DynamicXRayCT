function [V_true] = generateRigidVideo(N,nT,dt,options)
% Description: This function generates a rigid video with two balls 
% bouncing in a two-dimensional space.
%
% Inputs:
% N - Scalar, number of pixels along the x and y axes.
% nT - Scalar, number of time steps to simulate the bouncing of balls.
% dt - Scalar, time step size.
% options - Structure, containing parameter values for the balls.
%
% Outputs:
% I_true - 3D Matrix, the resulting video frames of the bouncing balls.
%
% Syntax:
% I_true = generateRigidVideo(N,nT,dt,options)
%
% Summary:
% The function takes four inputs, N, nT, dt, and options, and outputs a 
% video of the two balls bouncing in a two-dimensional space.


% Set default values for the inputs.
if nargin < 3, dt = 0.01; end
if nargin < 4, options = []; end

% Define the number of pixels along the x and y axes.
n = [N N];

% Get the parameters for the two balls from the options structure.
v1 = getoptions(options, 'v1', [0.5 1]);
v2 = getoptions(options, 'v2', [-0.8 0.5]);
c1 = getoptions(options, 'c1', [-0.5 -0.5]);
c2 = getoptions(options, 'c2', [-0.5 0.5]);

% Generate the spatial grid.
x = linspace(-1,1,n(1));
y = linspace(-1,1,n(2));
[xx,yy] = ndgrid(x,y);

% Initialize two balls.
I1 = zeros(n);
I2 = zeros(n);

% Initialize the first ball.
cx1 = c1(1);
cy1 = c1(2);
r1  = 0.2;
I1((xx-cx1).^2 + (yy-cy1).^2 <= r1^2) = 1;

% Initialize the second ball.
cx2 = c2(1);
cy2 = c2(2);
r2  = 0.15;
I2((xx-cx2).^2 + (yy-cy2).^2 <= r2^2) = 1;


% Initialize the velocity of the balls.
vx1 = v1(1);
vy1 = v1(2);
vx2 = v2(1);
vy2 = v2(2);

% Define the boundary of the two-dimensional space.
bDomx = zeros(n); bDomx([1:2 end-1:end],:) = 1;
bDomy = zeros(n); bDomy(:,[1:2 end-1:end]) = 1;

%% generate video

fprintf('\n'); 
fprintf('=============================================================\n');
fprintf('<strong>             Generating Phantom</strong>\n');
fprintf('                (Rigid Motion) \n');
fprintf('=============================================================\n');

V_true(:, :, 1) = I1 + I2;

% Define flags for balls hitting the boundary
I1xFlag = true;
I1yFlag = true;
I2xFlag = true;
I2yFlag = true;
I1I2Flag = true;

for t = 2:nT
    
    % Check if the balls hit the boundary and update velocity accordingly
    if vec(I1)' * vec(bDomx) <= 0
        I1xFlag = true;
    end
    if vec(I1)' * vec(bDomy) <= 0
        I1yFlag = true;
    end
    if vec(I2)' * vec(bDomx) <= 0
        I2xFlag = true;
    end
    if vec(I2)' * vec(bDomy) <= 0
        I2yFlag = true;
    end
    if vec(I1)' * vec(I2) <= 0
        I1I2Flag = true;
    end
    
    if vec(I1)' * vec(bDomx) > 0 && I1xFlag
        vx1 = -vx1;
        I1xFlag = false;
    end
    if vec(I1)' * vec(bDomy) > 0 && I1yFlag
        vy1 = -vy1;
        I1yFlag = false;
    end
    
    if vec(I2)' * vec(bDomx) > 0 && I2xFlag
        vx2 = -vx2;
        I2xFlag = false;
    end
    if vec(I2)' * vec(bDomy) > 0 && I2yFlag
        vy2 = -vy2;
        I2yFlag = false;
    end
    
    if vec(I1)' * vec(I2) > 0 && I1I2Flag
        vx1 = -vx1;
        vy1 = -vy1;
        vx2 = -vx2;
        vy2 = -vy2;
    end
    
    % Update ball 1
    cx1 = cx1 + vx1 * dt;
    cy1 = cy1 + vy1 * dt;
    I1 = zeros(size(I1));
    I1((xx - cx1).^2 + (yy - cy1).^2 <= r1^2) = 1;
    
    % Update ball 2
    cx2 = cx2 + vx2 * dt;
    cy2 = cy2 + vy2 * dt;
    I2 = zeros(size(I2));
    I2((xx - cx2).^2 + (yy - cy2).^2 <= r2^2) = 1;
    
    % Combine images
    V_true(:, :, t) = min(I1 + I2, 1);
end

% convert to single precision
V_true = single(V_true);

% Print information about the generated video
fprintf('Image size: %d x %d\n', [size(V_true, 1) size(V_true, 2)]);
fprintf('Timestamps: %d\n', nT);
fprintf('=============================================================\n');

end

