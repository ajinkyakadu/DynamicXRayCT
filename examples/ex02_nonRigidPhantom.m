%% Dynamic Tomography Script

% Clear variables and close all figures
clearvars; close all; clc;

%%% Set the seed for the random number generator
myStream = RandStream('mt19937ar', 'Seed', 1);
RandStream.setGlobalStream(myStream)

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
resDir = [pwd '/results/rigidPhantom/'];
if ~exist(resDir, 'dir'), mkdir(resDir); end

%% create a non-rigid-phantom
% converts the bell-phantom to tree-phantom using non-rigid deformations

N  = 512;       % Image dimensions
nT = 360;       % Number of frames,

% generate phantom
V_true = generateDeformationVideoA2B(N, nT, 'bell', 'tree');

%% setup the X-ray tomography (parallel-beam)

% define the angular step size and the noise level
dTheta = 5;
noiseL = 0.01;

% Perform the projections and store the measurements
[A, Af, y, theta] = singleShotProjections(V_true, dTheta, noiseL, use_cuda);

% get the size
n = size(V_true);

%% Algorithms

% Initialize parameters for a given phantom
maxIter= 1000;
params = getDefaultParameters(maxIter, dTheta, 'rigidPhantom');
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




