function [f, g, m, d] = levelSetDCT(x, funObj, use_cuda, options)
% LEVELSETDCT Computes the function, gradient, image, and boundary for a
%             level-set implicit representation of the binary images using 
%             the compressed discrete cosine basis. 
%
% DESCRIPTION:
%   This function implements the computation of the objective function and 
%   its gradient for the Level Set DCT problem. The problem is
%  
%       f(x) = F( m_1 h ( \Psi (x) ) )
% 
%   where \Psi denotes the discrete cosine transform, m_1 denotes the
%   intensity of the object, F denotes the objective function that relies 
%   on the image, and h is the Heaviside function (in our case, a smooth 
%   approximation of the Heaviside function). Here, the image is implicitly 
%   represented by the level-set of a function compressed by it's discrete 
%   cosine function (m = m_1 h(\Psi(x))). The gradient is obtained as:
% 
%       g(x) = m_1 \Psi^T d( \Psi (x) ) G( m_1 h ( \Psi (x) ) )
%
%   where G is the gradient of the objective function, \Psi^T is the
%   inverse discrete cosine transform, and d is the dirac-delta function
%   (in our case, an approximation of dirac-delta function).
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
%

% Set default values for options
if nargin < 4
    options = [];
end

% Retrieve options from the options structure
m1    = getoptions(options, 'm1', 1);
kappa = getoptions(options, 'kappa', 0.1);
n     = getoptions(options, 'n', []);
x_id  = getoptions(options, 'x_id', []);
sigma = getoptions(options, 'sigma', 0);

% Convert to the image domain
xv = zeros(prod(n), 1, 'single');
% if use_cuda
%     xv = gpuArray(xv);
% end
xv(x_id) = x;
xv = dctThree(reshape(xv, n), use_cuda);
xv = vec(xv);

% Heaviside width
hopt.epsi = kappa * (max(xv) - min(xv));

% Function and gradient
if sigma > 0
    [h, d, dp] = heavi(xv, hopt);
else
    [h, d] = heavi(xv, hopt);
end

% Scale the image by intensity
m = h * m1;

% Compute the objective and gradient wrt image
[f, g0] = funObj(m);
g0 = m1 * (g0 .* d);

% Compute boundary objective and corresponding gradient
if sigma > 0
    f = f + sigma * 0.5 * norm(d)^2;
    g0 = g0 + sigma * (dp .* d);
end

% Compute gradient wrt DCT coefficients
g0 = vec(idctThree(reshape(g0, n), use_cuda));
g  = g0(x_id);

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
