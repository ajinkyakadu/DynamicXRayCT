function [V] = generateDeformationVideoA2B(N, nT, phName1, phName2)
%GENERATEDEFORMATIONVIDEOA2B This function generates a deformed video given
%                           an input images and number of time steps.
%
% Summary: The generateDeformationVideo function generates a deformed video
% from an input image and number of time steps. The function takes two
% images, mainly for starting frame and end frame (both should be binary).
% It calculates the control points to deform starting image to final image,
% and then computes the deformation using 2nd order polynomial geometric
% transformation. The generates sequence of images are  stored in the 
% output variable V. The size of the video and the number of time steps are
% displayed in the console.
% 
% Inputs:
%   - N: input image size (default: 256)
%   - nT: The number of time steps (default: 360)
%   - phName1 - name of image for starting frame (default: 'bell')
%   - phName2 - name of image for last frame (default: 'tree')
%
% Outputs:
%   - V: The generated deformed video

% Set default options if not specified

if nargin < 1, N        = 256;    end
if nargin < 2, nT       = 360;    end
if nargin < 3, phName1  = 'bell'; end
if nargin < 4, phName2  = 'tree'; end

% Read in the input image 1
imageA = imread([phName1 '.png']);
imageA = single(imageA);
imageA = imresize(imageA, [N N],"nearest");
imageA = rescale(imageA);

% Read in the input image 2
imageB = imread([phName2 '.png']);
imageB = single(imageB);
imageB = imresize(imageB, [N N],"nearest");
imageB = rescale(imageB);


% Extract and match SIFT features between Image A and Image B
pointsA = detectSURFFeatures(imageA);
pointsB = detectSURFFeatures(imageB);

[featuresA, validPointsA] = extractFeatures(imageA, pointsA);
[featuresB, validPointsB] = extractFeatures(imageB, pointsB);

indexPairs     = matchFeatures(featuresA, featuresB, 'MatchThreshold', 100, ...
                    'MaxRatio', 0.99, 'Unique', true);
matchedPointsA = validPointsA(indexPairs(:, 1), :);
matchedPointsB = validPointsB(indexPairs(:, 2), :);

% Define the control points for the thin-plate spline transformation
controlPointsA = matchedPointsA.Location;
controlPointsB = matchedPointsB.Location;

fprintf('generated %d control points \n', size(controlPointsA,1));

%%% initialize
V(:,:,1)    = imageA;
warpedImage = imageA;

%%% loop over
for i = 2:nT

    % Calculate the control points for the current frame
    cA_fixed = controlPointsA + (i-1)/nT * (controlPointsB - controlPointsA);
    cA_moving = controlPointsA + i/nT * (controlPointsB - controlPointsA);
    
    % Calculate the transformation for the current frame
    tform_t = fitgeotform2d(cA_moving, cA_fixed, 'polynomial', 2);
    
    % Warp Image A towards Image B using the current transformation
    warpedImage = imwarp(warpedImage, tform_t, 'OutputView', imref2d(size(imageA)), 'FillValues', 0);
    
    % Assign the warped image to the current frame of the sequence
    V(:,:,i) = imbinarize(warpedImage);
end

% convert to single precision
V = single(V);

% print
fprintf('Image size: %d x %d \n',size(warpedImage));
fprintf('timestamps: %d \n',nT);

end