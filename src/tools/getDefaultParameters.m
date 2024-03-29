function params = getDefaultParameters(maxIter, dTheta, phantomName)
    % Initialize default parameters
    params = struct();
    
    % Common parameters
    params.common.numIter = maxIter;
    params.common.nB = floor(180 / dTheta);
    
    % Parameters for static binned reconstruction
    params.static.numIter = params.common.numIter;
    params.static.nB      = params.common.nB;

    % Parameters for Compressed Shape Sensing
    params.CSS.numIter = params.common.numIter;
    params.CSS.nB      = params.common.nB;

    % Parameters for Box + Temporal L2 Regularization
    params.boxL2.gamma    = 2.5e1; % Seems to be a good value for all sizes
    params.boxL2.numIter  = params.common.numIter;

    % Parameters for Dynamic Shape Sensing
    params.DSS.opt = struct('nBasis', [0.1, 0.1, 0.1], ... 
                               'xEstIter', 40, ...
                               'estimKp', 0, ...
                               'kpEstIter', 10, ...
                               'kappa', 0.1, ...
                               'tau', Inf);
    params.DSS.numIter       = floor(params.common.numIter/params.DSS.opt.xEstIter);
    params.DSS.plot_iterates = true;

    % Parameters for TV-TV-OF
    params.TVTVOF.opt = struct('alpha_TV', 0.0625 * 60, ...
                           'fac_alpha', 1, ...
                           'beta_TV', 1, ...
                           'fac_gamma', 1, ...
                           'numIterOuter', 10, ...
                           'numIterImage', 4000, ...
                           'numIterMotion', 1000);

    params.TVTVOF.motionVisuPara = struct('colorVisu', 'frame', ...
                                    'fps', 3);

    params.TVTVOF.ctInfo = struct('dim', 2, 'type', 'dynamic');
    
    % Adjust parameters based on the phantom name
    switch lower(phantomName)
        case 'cwidogtoy'
            params.DSS.opt.nBasis = [0.025, 0.025, 0.012];
        case 'rigidphantom'
            params.DSS.opt.nBasis = [0.05, 0.05, 0.05]; 
            params.DSS.opt.tau    = 1e4;
        case 'nonrigidphantom'
            params.DSS.opt.nBasis = [0.025, 0.025, 0.012];
    end
end
