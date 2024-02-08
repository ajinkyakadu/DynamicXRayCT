function X = dctThree(X, use_cuda)
% DCTThree - Perform Discrete Cosine Transform along three dimensions
%
% Syntax: X = dctThree(X)
%
% Inputs:
%   X - Input data to be transformed
%
% Outputs:
%   X - Transformed data
%
% The function performs the discrete cosine transform along three 
% dimensions using the specified 'Type' of transform (K = 2). The function 
% also casts the input to the host and back to the GPU to ensure 
% computation takes place on the CPU.

if nargin < 2, use_cuda = 0; end

% Type of transform
K = 1;

% Cast to host
if isa(X, 'gpuArray')
    X = gather(X);
end

% Perform DCT along each dimension
X = dct(X,[],1,'Type',K);
X = dct(X,[],2,'Type',K);
X = dct(X,[],3,'Type',K);

% Cast back to GPU
% if ~isa(X, 'gpuArray') && use_cuda
%     X = gpuArray(X);
% end


end


% Compute the 3D DCT using the 1D DCT function
function y = dct3(x)
    y = dct1(x);
    y = dct1(permute(y, [2 3 1]));
    y = dct1(permute(y, [2 3 1]));
    y = permute(y, [2 3 1]);
end

% Define a 1D DCT function for a gpuArray
function y = dct1(x)
    n = size(x,1);
    
    % Pad the input array with zeros
    % pad_size = 2^nextpow2(n+1) - n;
    % x_padded = [x; zeros(pad_size, size(x, 2), size(x,3))];
    
    % Compute the 1D DCT using the padded signal
    % y = fft(x_padded);
    y = fft(x, 2*n-1, 1);

    % Truncate the padded output to the original size of the input signal
    y = y(1:n, :, :);

    % Scale the DCT by a factor of sqrt(2/n)
    y = 2*real(y)/sqrt(2*n-1);

    y(1,:,:) = y(1,:,:) / sqrt(2);
end