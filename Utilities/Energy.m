function [ElasticEnergy,BendingEnergy] = Energy(BorderCoords)
% -------------------------------------------------------------------------
% Function that calculates the elastic and bending enery of coordinates of
% a point cloud. These coordinates should be border points of the point 
% cloud of the calculation does not have much meaning. These coordinates
% are also sequential coordinates that follow the trajectory of the border.
% Example on how to use it:
%   [EE,BE] = Energy(PointCloudBorders);
% -------------------------------------------------------------------------
% Input:
%   BorderCoords:   The coordinates of the border (must be ordered).
%                   This should be a 2-column matrix (x in the first column
%                   and y in the second colum).
%
% Output:
%   ElasticEnergy:  The elastic energy of the border.
%   BendingEnergy:  A cell containing the polyshapes of the individual
%                   holes.
% -------------------------------------------------------------------------
% Code written by:
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
%   Siewert Hugelier    Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   siewert.hugelier@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   Hugelier, S., Kim, H., Gyparaki, M.T., Bond, C., Tang, Q., 
%   Santiago-Ruiz, A.N., Porta, S., Lakadamyali, M. ECLiPSE: a versatile 
%   classification technique for structural and morphological analysis of 
%   super-resolution microscopy data. BioRxiv (2023). 
%   DOI: https://doi.org/10.1101/2023.05.10.540077.
% -------------------------------------------------------------------------

% Check if the input data is in the correct format. If not, show an error.
if ~ismatrix(BorderCoords) || size(BorderCoords,2) ~= 2 || size(BorderCoords,1) <= 2
    error("The input should be specified as a 2-column matrix (x in colum 1, and y in column 2) that contains more than 2 input coordinates.");
end

% Calculate the first and second derivative of the coordinates. Do this in 
% a columwise fashion.
Diff1 = diff(BorderCoords,1,1);
Diff2 = diff(BorderCoords,2,1);

% Calculate the elastic and bending energy.
ElasticEnergy = 0.5 * sum(Diff1(:).^2); % Elastic energy (1/2 * Δx²).
BendingEnergy = 0.5 * sum(Diff2(:).^2); % Bending energy (1/2 * ΔΔx²).

end

