function results = computeAlgorithmMetrics(runFlags, algoResults, V_true)
% runFlags: Struct with fields indicating whether each algorithm should be run (e.g., runFlags.static)
% algoResults: Struct containing the results and metadata for each algorithm (e.g., algoResults.static.V, algoResults.static.time)
% V_true: The ground truth image for comparison

% Initialize results structure
results = struct();
algorithms = {'static', 'CSS', 'boxL2', 'DSS', 'TVTVOF'};

fprintf('\n'); 
fprintf('=============================================================\n');
fprintf('<strong>                  Algorithm Metrics </strong>\n');
fprintf('                (PSNR, SSIM, Dice, etc) \n');
fprintf('=============================================================\n');


for algo = algorithms
    algoKey = algo{1};
    if isfield(runFlags, algoKey) && runFlags.(algoKey)
        % Check if the algorithm results exist
        if isfield(algoResults, algoKey) && isfield(algoResults.(algoKey), 'V')
            V_algo = algoResults.(algoKey).V;
            t_algo = algoResults.(algoKey).time;
            energy_algo = isfield(algoResults.(algoKey), 'energy') * algoResults.(algoKey).energy;

            % Compute metrics
            psnr_val = psnr(V_algo, V_true);
            ssim_val = ssim(V_algo, V_true);
            dice_val = diceCoefficient(V_algo, V_true);

            % Store results
            results.(algoKey).rec    = V_algo;
            results.(algoKey).psnr   = psnr_val;
            results.(algoKey).ssim   = ssim_val;
            results.(algoKey).dice   = dice_val;
            results.(algoKey).time   = t_algo;
            results.(algoKey).energy = energy_algo;

            % Print results
            fprintf('%6s: PSNR:%.4f | SSIM:%.4f | Dice:%.4f | time: %.4f \n', algoKey, ...
                psnr_val, ssim_val, dice_val, t_algo);
        else
            fprintf('Results for %s not available.\n', algoKey);
        end
    end
end

fprintf('=============================================================\n');

end
