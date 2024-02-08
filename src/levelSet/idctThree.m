function X = idctThree(X, use_cuda)
% IDCTThree - Perform Inverse Discrete Cosine Transform along three 
% dimensions
%
% Syntax: X = idctThree(X)
%
% Inputs:
%   X - Input data to be transformed
%
% Outputs:
%   X - Transformed data
%
% The function performs the inverse discrete cosine transform along three 
% dimensions using the specified 'Type' of transform (K = 2). The function 
% also casts the input to the host and back to the GPU to ensure 
% computation takes place on the CPU.

if nargin < 2, use_cuda = 0; end

% Type of transform
K = 1;

% Cast to host
if use_cuda
    X = gather(X);
end

% % Perform inverse DCT along each dimension
X = idct(X,[],3,'Type',K);
X = idct(X,[],2,'Type',K);
X = idct(X,[],1,'Type',K);

% % Cast back to GPU
% if use_cuda
%     X = gpuArray(X);
% end

% X = idct3(X);

end


% Compute the 3D IDCT using the 1D IDCT function
function y = idct3(x)
    y = idct1(x);
    y = idct1(permute(y, [2 3 1]));
    y = idct1(permute(y, [2 3 1]));
    y = permute(y, [2 3 1]);
end

% Define a 1D IDCT function for a gpuArray
function y = idct1(x)
    n = size(x,1);

    % Scale the input by a factor of sqrt(2/n)
    x(1,:,:) = x(1,:,:) / sqrt(2);

    % Compute the 1D IDCT using the padded signal
    y_padded = ifft(x, 2*n-1, 1);

    % Truncate the padded output to the original size of the input signal
    y = y_padded(1:n, :, :);

    % Scale the IDCT by a factor of sqrt(2/n)
    y = 2*real(y)*sqrt(2*n-1); % 
end

% % Define a 1D IDCT function for a gpuArray
% function y = idct1(x)
%     n = size(x,1);
%     x(1,:,:) = x(1,:,:) / sqrt(2);
%     x = [x; conj(x(end-1:-1:2,:,:))];
%     y = ifft(x);
%     y = real(y(1:n,:,:));
%     y = y * sqrt(2/n);
% end