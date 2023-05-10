function Curvature = Curvature(Points)
% -------------------------------------------------------------------------
% Function that calculates the curvature of points. These points are most
% likely border points.
% Example on how to use it:
%   CurvaturePointCloud1 = Curvature(BorderPoints_PointCloud1);
% -------------------------------------------------------------------------
% Input:
%   Points: The (border)points of the curve of which the curvature has to
%           be calculated. This should be a (n x 2) or (n x 3) vector. If
%           it is an (n x 2), then the 3rd dimension will be set at 0.
%
% Output:
%   Curvature:  The curvature of each point (determined by the radius of
%               the oscillating circle using the point and its neighbouring
%               points.
% -------------------------------------------------------------------------
% Code written by:
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
[M,N] = size(Points); % Determine the size of the coordinates.
if N ~= 2 && N ~= 3
    error("The input should be specified as a 2-column or 3-column matrix (x in colum 1, y in column 2, and z in column 3)");
end

% Add a third column for the z-position to the coordinates if they are 2D.
if N == 2
    Points(:,3) = 0;
end

% Pre-allocate vectors for memory purposes.
Radius = cell(M,1); Radius(:) = {NaN}; % Pre-allocate for speed reasons.

% Set up the vector containing the neighbouring point to make up the
% triangles to calculate the radius of the curvature. Check also if the
% curve is closed, and add those neighbouring points manually. Convert the triangles 
% to cells to make calculations faster, and then extract the coordinates
% associated to these indices.
Triangles = [(1:M-2)' (2:M-1)' (3:M)']; % Make the vector of neighbouring points: {i-1; i; i+1}.

if norm(Points(1,:)-Points(M,:)) < 1E-8 % Check whether or not the curve is closed.
    Triangles = vertcat([M-1 1 2],Triangles,[M-1 M 2]); % If the curve is close, then add the start and end point as well.
    Triangles = mat2cell(Triangles,ones(M,1),3); % Convert to a cell.
else
    Triangles = mat2cell(Triangles,ones(M-2,1),3); % Convert to a cell.
end

CoordsTriangles = cellfun(@(x) Points(x,:)',Triangles,'UniformOutput',false); % Extract the coordinates of the points for each triangle.

% Calculate the radius of the curvature.
Crossproduct = cellfun(@(x) norm(cross(x(:,1)-x(:,2),x(:,2)-x(:,3))),CoordsTriangles,'UniformOutput',false); % Magnitude of the crossproduct between the two vectors describing the direction.
NormLeft = cellfun(@(x) norm(x(:,2)-x(:,1)),CoordsTriangles,'UniformOutput',false); % Euclidean distance between the coordinate of interest and the previous coordinate.
NormRight = cellfun(@(x) norm(x(:,2)-x(:,3)),CoordsTriangles,'UniformOutput',false); % Euclidean distance between the coordinate of interest and the next coordinate.
NormLeftRight = cellfun(@(x) norm(x(:,1)-x(:,3)),CoordsTriangles,'UniformOutput',false); % Euclidean distance between the previous coordinate and the next coordinate.

if norm(Points(1,:)-Points(M,:)) < 1E-8
    Radius = cellfun(@(w,x,y,z) {w*x*y/(2*z)},NormLeftRight,NormLeft,NormRight,Crossproduct); % Radius of the curvature (when curve is closed).
else
    Radius(2:M-1) = cellfun(@(w,x,y,z) {w*x*y/(2*z)},NormLeftRight,NormLeft,NormRight,Crossproduct); % Radius of the curvature (when curve is open).
end
Radius = vertcat(Radius{:}); % Vectorize the cell.

% Calculate the curvature
Curvature = 1 ./ Radius; % curvature is 1 over the radius.

end