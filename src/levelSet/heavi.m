function [h, d, dp] = heavi(x, options)

% Check if options is provided, if not, set options to an empty structure
if nargin < 2
    options = struct();
end

% Get options from input
type = getoptions(options, 'type', 'compact');
epsi = getoptions(options, 'epsi', 1e-16);
thr  = getoptions(options, 'thr', 0);

% Shift x by threshold
x = x - thr;

% Switch calculation based on type of heaviside function
switch type
    case 'global'
        h = 0.5 * (1 + 2 / pi * atan(pi * x / epsi));
        
        if nargout > 1
            divisor = (epsi * ((x.^2 * pi^2) / epsi^2 + 1));
            d = 1 ./ divisor;
        end
        
        if nargout > 2
            dp = - (2 * pi^2 / epsi) * (x ./ (divisor.^2));
        end
        
    case 'compact'
        h = 0*x;
        id = find(abs(x) < epsi);
        
        h(id) = 0.5 * (1 + x(id) / epsi + (1/pi) * sin(pi * x(id) / epsi));
        h(x >= epsi)  = 1;
        h(x <= -epsi) = 0;
        
        if nargout > 1
            d = 0*x;
            d(id) = 0.5 * (1 / epsi) * (1 + cos(pi * x(id) / epsi));
        end
        
        if nargout > 2
            dp = 0*x;
            dp(id) = - (pi / (2 * epsi^2)) * sin(pi * x(id) / epsi);
        end
end

end