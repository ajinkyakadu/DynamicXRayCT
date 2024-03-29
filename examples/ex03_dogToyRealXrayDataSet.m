%% Dynamic Tomography Script

% Clear variables and close all figures
clearvars; close all; clc;


% check whether we are in display mode
visulization = usejava('desktop');

%%% flags
% Algorithm execution flags
runFlags = [];
runFlags.static = false;
runFlags.CSS    = false;
runFlags.boxL2  = true;
runFlags.DSS    = true;
runFlags.TVTVOF = false;

% flag to use CUDA compute (read memory around 12 GB)
use_cuda = true;

if use_cuda
    reset(gpuDevice());
end

%%% save directory
resDir = [pwd '/results/CWIdogToyCT/'];
if ~exist(resDir, 'dir'), mkdir(resDir); end

%% (download and) load dataset
% Dataset is available on Zenodo (This is the processed dataset; hence can
% be included and work with directly).


% URL of the dataset
url = 'https://zenodo.org/records/10808364/files/GrayBone90kV4FilterPreprocessed.mat?download=1';

% Specify the file name
fileName = 'GrayBone90kV4FilterPreprocessed.mat';
filePath = fullfile(resDir, fileName); % Full path to the file

% Check if the file already exists; If the file does not exist, download 
% the file
if exist(filePath, 'file')
    disp(['File already exists at: ', filePath]);
else
    
    websave(filePath, url);
    disp(['File downloaded to: ', filePath]);
end


% Load the dataset into MATLAB workspace
load(filePath);
disp('Dataset loaded into MATLAB workspace here');


%% setup the X-ray tomography (parallel-beam)

% Define the angular step size
nT     = 600;
dTheta = 5;

[A, Af, y, angles] = getRealTomographicMeasurements(stack, nT, dTheta, use_cuda);

V_true = im2single(stack.rec(:,:,1:nT));

clear stack


%% Algorithms

% Initialize parameters for a given phantom
maxIter= 1000;
params = getDefaultParameters(maxIter, dTheta, 'cwidogtoy');
algoResults = struct();

%%% Static binned reconstruction
if runFlags.static 
    
    [V_static, en_static, t_static] = computeStaticSolution(Af, y, n, ...
        params.static.nB, params.static.numIter, use_cuda, false);

    algoResults.static.V      = V_static;
    algoResults.static.time   = t_static;
    algoResults.static.energy = en_static;
end

%%% Compressed Shape Sensing
if runFlags.CSS

    [V_CSS, en_CSS, t_CSS] = computeCompressedShapeSensingSolution(Af, y, ...
        n, params.CSS.nB, params.CSS.numIter, use_cuda);

    algoResults.CSS.V      = V_CSS;
    algoResults.CSS.time   = t_CSS;
    algoResults.CSS.energy = en_CSS;
end

%%% Box + Temporal L2 Regularization
if runFlags.boxL2

    [V_boxL2, en_boxL2, t_boxL2] = computeBoxL2RegSolution(A, y(:), n, params.boxL2.gamma, ...
        params.boxL2.numIter, use_cuda);

    algoResults.boxL2.V      = V_boxL2;
    algoResults.boxL2.time   = t_boxL2;
    algoResults.boxL2.energy = en_boxL2;
end

%%% Dynamic Shape Sensing
if runFlags.DSS

    [V_DSS, en_DSS, t_DSS, f_DSS] = computeDynamicShapeSensingSolution(A, y(:), n, ...
        params.DSS.numIter, use_cuda, params.DSS.opt, V_boxL2, params.DSS.plot_iterates);

    algoResults.DSS.V      = V_DSS;
    algoResults.DSS.time   = t_DSS;
    algoResults.DSS.energy = en_DSS;
end


%%% TV-TV-OF
if runFlags.TVTVOF
    
    [V_TVTVOF, U_TVTVOF, t_TVTVOF, info_TVTVOF] = ...
        computeTVTVOFSolution(Af, y, n, params.TVTVOF.fac_alpha * params.TVTVOF.alpha_TV, ...
        params.TVTVOF.beta_TV, params.TVTVOF.fac_gamma * params.TVTVOF.gamma_L2, ...
        params.TVTVOF.numIterOuter, params.TVTVOF.numIterImage, params.TVTVOF.numIterMotion);

    algoResults.TVTVOF.V      = V_TVTVOF;
    algoResults.TVTVOF.time   = t_TVTVOF;
    algoResults.TVTVOF.energy = en_TVTVOF;
end

%% Compute Algorithm Metrics

results = computeAlgorithmMetrics(runFlags, algoResults, V_true);

%% Save Results

execTime = char(datetime("now", "Format","uuuuMMdd'T'HHmmss"));
save([resDir 'results_ ' execTime '.mat'], 'params', 'results', '-v7.3');




