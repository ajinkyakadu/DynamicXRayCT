function [diceCoefficient] = diceCoefficient(reconArray, groundTruthArray)
% Calculates Dice coefficient between two binary images
%
% Inputs:
%   reconArray : reconstructed binary array
%   groundTruthArray : ground truth binary array
%
% Output:
%   diceCoefficient : Dice coefficient between two binary arrays

% Convert input arrays into column vectors
reconArray = vec(reconArray);
groundTruthArray = vec(groundTruthArray);

% If reconstructed array has more than two unique values, apply Otsu's
% thresholding method to convert it into binary array
if length(unique(reconArray)) > 2
t = multithresh(reconArray);
reconArray = reconArray > t;
end

% Calculate the intersection, sum of reconstructed array and sum of ground
% truth array
intersection = sum(reconArray & groundTruthArray);
reconstructed_sum = sum(reconArray);
ground_truth_sum = sum(groundTruthArray);

% Calculate the Dice coefficient
diceCoefficient = (2 * intersection) / (reconstructed_sum + ground_truth_sum);

end