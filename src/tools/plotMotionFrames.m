function [] = plotMotionFrames(img1, img2)

n = size(img1);
blk_size = [7 7];

hbm = vision.BlockMatcher('ReferenceFrameSource',...
    'Input port','BlockSize',blk_size);
hbm.OutputValue = 'Horizontal and vertical components in complex form';
halphablend = vision.AlphaBlender;

motion = hbm(img1,img2);

img12 = halphablend(img2,img1);

% Compute the gradient of the image
[Gx, Gy] = imgradientxy(img1);
grad_mag = sqrt(Gx.^2 + Gy.^2);

% Threshold the gradient image to obtain a binary edge map
edge_map = grad_mag > 0.1;

% Dilate the edge map to expand the boundary by a few pixels
se = strel('disk', 3);
edge_map_dilated = imdilate(edge_map, se);

% Use the dilated edge map to mask the motion vectors
[X,Y] = meshgrid(1:blk_size(1):size(img1,2),1:blk_size(2):size(img2,1));

edge_map_dilated = imresize(edge_map_dilated, size(X), "nearest");

masked_motion_x = real(motion) .* edge_map_dilated;
masked_motion_y = imag(motion) .* edge_map_dilated;

imshow(img12);
hold on

% Calculate the angle of each vector in radians
angles = atan2(masked_motion_y, masked_motion_x);

% Normalize the angles to [0, 1] range
normalized_angles = (angles + pi) / (2 * pi);

% Create an HSV colormap
cmap = hsv(256);

% Map the normalized angles to the corresponding colormap index
color_indices = round(normalized_angles * (size(cmap, 1) - 1)) + 1;

% Assign colors to the vectors based on the colormap
vector_colors = cmap(color_indices, :);

% Plot the vector field with colors from the colormap
idxSelectX = find(masked_motion_x);
idxSelectY = find(masked_motion_y);

idxSelect = [idxSelectX;idxSelectY];
idxSelect = unique(idxSelect);

AutoScaleFactor = 1.5;
for i = 1:length(idxSelect)
    quiver(X(idxSelect(i)), Y(idxSelect(i)), masked_motion_x(idxSelect(i)),...
        masked_motion_y(idxSelect(i)), 0, 'Color', vector_colors(idxSelect(i),:), ...
        'MaxHeadSize', 1, 'AutoScale', 'off', 'LineWidth', 1);
end

hold off

end
