function [fig] = computeTSNEplot(V)

n = size(V);
V = reshape(V, [n(1)*n(2), n(3)]);

options = struct();
options.MaxIter = 10000;
options.TolFun  = 1e-10;

rng default; % For reproducibility
reduced_data = tsne(V', 'NumDimensions', 2, ...
    'Distance', 'correlation', 'Algorithm','exact', 'LearnRate',500, ...
    'Options',options, 'Standardize', false, 'Perplexity', 50,...
    'Exaggeration', 10);


% Modify the jet colormap to change brightness
n_colors = n(3);
cmap     = jet(n_colors);

fig = figure;
hold on;

% Interpolate the data points using a smooth curve
t = 1:n(3);
nP = 1000;
tq = linspace(1, n(3), nP);

smooth_x = interp1(t, reduced_data(:, 1), tq, 'pchip');
smooth_y = interp1(t, reduced_data(:, 2), tq, 'pchip');
color_vals = linspace(1, n(3), nP);

% Plot smooth curve
surface([smooth_x(:), smooth_x(:)], [smooth_y(:), smooth_y(:)], ...
    [zeros(size(color_vals(:))), zeros(size(color_vals(:)))], ...
    [color_vals(:), color_vals(:)], ...
    'EdgeColor', 'interp', 'Marker', 'none', 'LineWidth', 2);

% Define text offset
offset = 0.05;

% Calculate axis range
yrange = max(reduced_data(:, 2)) - min(reduced_data(:, 2));
xrange = max(reduced_data(:, 1)) - min(reduced_data(:, 1));


% Add text labels for starting, middle, and end points
text(reduced_data(1, 1) + offset*xrange, reduced_data(1, 2) + offset*yrange, 't = 1', 'FontSize', 12, 'FontWeight', 'bold');
text(reduced_data(n(3)/2, 1) + offset*xrange, reduced_data(n(3)/2, 2) + offset*yrange, sprintf('t = %d', n(3)/2), 'FontSize', 12, 'FontWeight', 'bold');
text(reduced_data(end, 1) + offset*xrange, reduced_data(end, 2) + offset*yrange, sprintf('t = %d', n(3)), 'FontSize', 12, 'FontWeight', 'bold');


% Mark starting, middle, and end points on the plot
scatter(reduced_data(1, 1), reduced_data(1, 2), 60, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b', 'LineWidth', 1.5);
scatter(reduced_data(n(3)/2, 1), reduced_data(n(3)/2, 2), 60, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'LineWidth', 1.5);
scatter(reduced_data(end, 1), reduced_data(end, 2), 60, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'LineWidth', 1.5);

% Set axis limits with margin
x_margin = 0.3 * range(reduced_data(:, 1));
y_margin = 0.3 * range(reduced_data(:, 2));
xlim([min(reduced_data(:, 1)) - x_margin, max(reduced_data(:, 1)) + x_margin]);
ylim([min(reduced_data(:, 2)) - y_margin, max(reduced_data(:, 2)) + y_margin]);

% title('Low-dimensional Temporal Manifold');
grid on; box on; axis square;
colormap(cmap);
% caxis([1, n(3)]);
set(gca,'xtick',[]);
set(gca,'ytick',[]);

% Determine the axis limits
x_min = min(reduced_data(:, 1)) - x_margin;
x_max = max(reduced_data(:, 1)) + x_margin;
y_min = min(reduced_data(:, 2)) - y_margin;
y_max = max(reduced_data(:, 2)) + y_margin;

% Add arrows for the axes
arrow_length_x = 0.1 * (x_max - x_min);
arrow_length_y = 0.1 * (y_max - y_min);
offset = 0.05;
quiver(x_min+offset*xrange, y_min+offset*yrange, arrow_length_x, 0, 'k', 'LineWidth', 1, 'MaxHeadSize', 0.5);
quiver(x_min+offset*xrange, y_min+offset*yrange, 0, arrow_length_y, 'k', 'LineWidth', 1, 'MaxHeadSize', 0.5);

% Add text labels for the axes
text_offset = 0.01;
text(x_min + offset*xrange + arrow_length_x, y_min + offset*yrange + (y_max - y_min) * text_offset, 't-SNE1', 'FontSize', 12);
text(x_min + offset*xrange + (x_max - x_min) * text_offset, y_min + offset*yrange + arrow_length_y, 't-SNE2', 'FontSize', 12);

hold off;
pause(.001);

end