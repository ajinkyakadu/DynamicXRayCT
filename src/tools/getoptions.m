function v = getoptions(options, name, v, mandatory)
% Get Options Function - The function retrieves an option parameter from
% the input options structure.
%
% Inputs:
%   options - A structure containing the options and their values
%   name - The name of the option to retrieve
%   v - The default value of the option
%   mandatory (optional) - A flag indicating whether the option is 
%               mandatory or not (default: 0)
%
% Outputs:
%   v - The retrieved value of the option

% Check if the mandatory flag was provided
if nargin < 4
    mandatory = 0;
end

% Check if the option exists in the options structure
if isfield(options, name)
    v = eval(['options.' name ';']);
elseif mandatory
    % Raise an error if the option is mandatory and not provided
    error(['You have to provide options.' name '.']);
end

end