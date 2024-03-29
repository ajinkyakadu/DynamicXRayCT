function [] = gatherTSNEplots(runFlags, algoResults, V_true, figDir)
% runFlags: Struct with fields indicating whether each algorithm should be run (e.g., runFlags.static)
% algoResults: Struct containing the results and metadata for each algorithm (e.g., algoResults.static.V, algoResults.static.time)
% V_true: The ground truth image for comparison

algorithms = {'static', 'CSS', 'boxL2', 'DSS', 'TVTVOF'};

fprintf('\n'); 
fprintf('=============================================================\n');
fprintf('<strong>                  Computing TSNE plots </strong>\n');
fprintf('=============================================================\n');


% Plot slices of V_true for comparison
[fig1] = computeTSNEplot(V_true);
export_fig(fullfile(figDir, 'GT_tsne'), '-png', '-m1', fig1);
close(fig1)

% Counter for subplot positioning
for algo = algorithms
    algoKey = algo{1};
    if isfield(runFlags, algoKey) && runFlags.(algoKey)
        % Check if the algorithm results exist
        if isfield(algoResults, algoKey) && isfield(algoResults.(algoKey), 'V')
            fprintf('Running for %s ...\n', algoKey);
            V_algo = algoResults.(algoKey).V;
            [fig1] = computeTSNEplot(V_algo);
            export_fig(fullfile(figDir, [algoKey '_tsne']), '-png', '-m1', fig1);
            close(fig1)
        else
            fprintf('Results for %s not available.\n', algoKey);
        end
    end
end

fprintf('=============================================================\n');

end


