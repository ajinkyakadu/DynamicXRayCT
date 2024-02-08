function setup()

% Define toolboxes and download links
toolboxes = {
    'ASTRA', 'https://github.com/astra-toolbox/astra-toolbox/releases/download/v2.1.0/astra-toolbox-2.1.0-matlab-win-x64.zip', 'astra*';
    'SPOT', 'https://www.cs.ubc.ca/labs/scl/spot/spot-1.2.zip', 'spot*';
    'FelixMatlabTools', 'https://github.com/FelixLucka/FelixMatlabTools/archive/refs/tags/v24.02.01.zip', 'FelixMatlabTools*';
    'FelixMatlabSolver', 'https://github.com/FelixLucka/FelixMatlabSolvers/archive/refs/tags/v24-02-01.zip', 'FelixMatlabSolvers*';
    'MinConf', 'https://www.cs.ubc.ca/~schmidtm/Software/minConf.zip', 'minConf*';
    'iFEM', 'https://github.com/lyc102/ifem/archive/refs/heads/master.zip', 'iFEM*';
};

% Create target folder for toolboxes
home_folder = pwd();

target_folder = fullfile(home_folder, 'src', 'external');
if ~exist(target_folder, 'dir')
    mkdir(target_folder);
end

fprintf('=============================================================\n');
fprintf('               setting up DynamicXRayCT                      \n');
fprintf('-------------------------------------------------------------\n');

% Download and extract toolboxes
for i = 1:size(toolboxes, 1)
    toolbox_name = toolboxes{i, 1};
    toolbox_url = toolboxes{i, 2};
    toolbox_pattern = toolboxes{i, 3}; % Pattern to check existence

    % Get a list of directories before extraction
    pre_extraction_dirs = dir(fullfile(target_folder, toolbox_pattern));
    
    
    % Check if the toolbox is already installed by searching for a matching pattern
    dir_contents = dir(fullfile(target_folder, toolbox_pattern));
    if isempty(dir_contents)
        fprintf('Downloading %s...\n', toolbox_name);
        
        % Download the ZIP file
        zip_filename = fullfile(target_folder, [toolbox_name, '.zip']);
        websave(zip_filename, toolbox_url);
        
        % Extract the ZIP file
        fprintf('Extracting %s...\n', toolbox_name);
        unzip(zip_filename, target_folder);

        % Delete the ZIP file
        delete(zip_filename);

        % Get a list of directories after extraction
        post_extraction_dirs = dir(fullfile(target_folder, toolbox_pattern));
        
        % Find the new directory (or directories) added by the extraction
        new_dirs = setdiff({post_extraction_dirs.name}, {pre_extraction_dirs.name});
        
        extracted_toolbox_path = fullfile(target_folder, new_dirs{1});
        
        fprintf('%s successfully installed.\n', toolbox_name);
    else
        fprintf('%s is already installed.\n', toolbox_name);
    end
end

% Add the target folder and its subfolders to the MATLAB path
addpath(genpath(fullfile(home_folder, 'src')));

fprintf('All the toolboxes and src files added to Path.\n');
fprintf('Toolbox setup complete.\n');

fprintf('=============================================================\n');

end
