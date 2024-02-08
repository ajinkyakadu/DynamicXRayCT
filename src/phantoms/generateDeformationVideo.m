function [V] = generateDeformationVideo(N, nT, phName, options)
%GENERATEDEFORMATIONVIDEO This function generates a deformed video given an
%                         input image and number of time steps.
%
% Summary: The generateDeformationVideo function generates a deformed video
% from an input image and number of time steps. The function has options 
% for the maximum rotation, translation, and scale to be applied to the 
% video. The function generates random rotation angles, translation 
% vectors, and scaling values. It then applies rotation, translation, and 
% binarization to each frame in the video. The final deformed video is 
% stored in the output variable V. The size of the video and the number of 
% time steps are displayed in the console.
% 
% Inputs:
%   - N: input image size 
%   - nT: The number of time steps
%   - options: A struct of options (optional)
%       rot_max: The maximum rotation (default: 1)
%       trans_max: The maximum translation (default: 1)
%       scale_max: The maximum scale (default: 0.01)
%
% Outputs:
%   - V: The generated deformed video

% Set default options if not specified

if nargin < 3, phName  = 'foam'; end
if nargin < 4, options = []; end

% get parameters
freq         = getoptions(options, 'freq',      10);
amp          = getoptions(options, 'amp',       2);
close_rad    = getoptions(options, 'close_rad', 10);
erode_rad    = getoptions(options, 'erode_rad', 5);
applyMorphOp = getoptions(options, 'apply_morphOp', false);


% Read in the input image
% Read the image and resize it
inputImage = imread([pwd '/images/' phName '.png']);
inputImage = single(inputImage);
inputImage = imresize(inputImage, [N N]);
inputImage = rescale(inputImage);

% Define the points for the initial and final meshes
[X,Y] = meshgrid(1:size(inputImage,2),1:size(inputImage,1));

% Define the Lucas-Kanade parameters
windowSize       = 15;
numPyramidLevels = 3;

V = zeros([N N nT], 'single');

Xn = X;
Yn = Y;


% Loop through the frames and generate the non-rigid deformation
for i = 1:nT

    % Calculate amplitude for current frame
    amp_t = (1 - 0.5*(i-1)/(nT-1)) * amp;

    % Warp the input image using the current mesh
    Xn = Xn + amp_t*cos(Y/freq) + amp_t/2*sin(2*X/freq).*cos(2*Y/freq);
    Yn = Yn + amp_t*sin(X/freq) + amp_t/2*sin(2*Y/freq).*cos(2*X/freq);
    
    warpedImage = interp2(X,Y,inputImage,Xn,Yn,'nearest',0);
    
    % binarize image
    warpedImage = imbinarize(warpedImage);
    
    if applyMorphOp
        % close image
        if close_rad > 0
            se_close = strel('disk', close_rad);
            warpedImage = imclose(warpedImage, se_close);
        end

        % erode small contents
        if erode_rad > 0
            se_erode = strel('disk', erode_rad);
            warpedImage = imerode(warpedImage, se_erode);
        end
    end
   
    % assign
    V(:,:,i) = warpedImage;
end

fprintf('Image size: %d x %d \n',size(warpedImage));
fprintf('timestamps: %d \n',nT);
    
end
