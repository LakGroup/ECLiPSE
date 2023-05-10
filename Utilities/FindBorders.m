function [BorderCoords, PolyshapeCloud] = FindBorders(Cloud)
% -------------------------------------------------------------------------
% Function that finds the inner and outer borders of a point cloud. The
% inner coordinates are only there when 'holes' are present in the data.
% The function uses a combination of the alphaShape and polyshape functions
% of Matlab and then extract the borders sequentially (they are ordered
% according to their trajectory).
% Example on how to use it:
%   PointCloud_AS = alphaShape(PointCloudCoords);
%   [PointCloud_BorderCoordinates, PointCloud_Polyshape] = FindBorders(PointCloud_AS);
% -------------------------------------------------------------------------
% Input:
%   Cloud:  The point cloud of the data.
%           This can either be a 2-column matrix (x in the first column
%           and y in the second colum), or an alphaShape (see 'help
%           alphaShape').
%
% Output:
%   BorderCoords:   A cell containing the coordinates of the outer border
%                   and the coordinates of each of the holes in the cloud 
%                   (if any) (i.e., the inner borders).
%   ShapePolygon:   A polyshape of the combined polygon (a direct
%                   'translation' of the alphaShape).
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

% Turn warnings off, to not flood the command window
warning('off','all');

% Check if the input data is in the correct format. If not, show an error.
if ~(ismatrix(Cloud) || isa(Cloud,'alphaShape')) || size(Cloud,2) > 2
    error("The input should be specified as a 2-column matrix (x in colum 1, and y in column 2) or as an alphaShape (see 'help alphaShape')");
end

% Check if the input is a matrix or an alphaShape, and calculate the
% alphaShape if needed.
if isa(Cloud,'double')
    Cloud = alphaShape(Cloud); % Create the alphaShape (uses Delaunay triangulation).
    AlphaOneRegion = criticalAlpha(Cloud,'one-region'); % Calculate the alpha corresponding to a 'one region' alphaShape, to make sure that all (x,y) coordinates are connected.
    Cloud.Alpha = AlphaOneRegion; % Apply it to the alphaShape.
end

% Get the outer boundary of the alphaShape and make a polyshape of it.
I = boundary(Cloud.Points,1); % Calculate the boundary points of the alphaShape and retrieve their indices.
BorderPoints = Cloud.Points(I,:); % Retrieve the actual boundary points.
PolygonCluster = polyshape(BorderPoints(:,1),BorderPoints(:,2)); % Create a polygon of the boundary points.

% Make the filled alphaShape of the original alphaShape.
Cloud_Filled = Cloud; % Copy the original alphaShape (to keep the same Alpha value).
Cloud_Filled.HoleThreshold = 1E15; % Set the hole threshold very high so all holes are filled.

% Extract the Delauny triangles from the two different alphaShape clouds.
% Then compare the coordinates of these triangles to see which ones belong 
% to the holes (the indices cannot be used here because they change between
% the two different alphaShapes).
[IndivTriangles,IndivPoints] = alphaTriangulation(Cloud); % Extract the different Delauny triangles of the alphaShape with holes.
[IndivTriangles_Filled,IndivPoints_Filled] = alphaTriangulation(Cloud_Filled); % Extract the different Delauny triangles of the alphaShape without holes.

IndivTriangles = num2cell(IndivTriangles); % Convert the triangle indices to a cell for each triangle (easier to work with).
IndivTriangles_Filled = num2cell(IndivTriangles_Filled); % Convert the triangle indices to a cell for each triangle (easier to work with).

TriangleCoordinates = cellfun(@(x) IndivPoints(x,:),IndivTriangles,'UniformOutput',false); % Extract the coordinates of each Delauny triangle of the alphaShape with holes.
TriangleCoordinates = cell2mat(TriangleCoordinates); % Convert to a matrix (1 row per triangle).
TriangleCoordinates_Filled = cellfun(@(x) IndivPoints_Filled(x,:),IndivTriangles_Filled,'UniformOutput',false); % Extract the coordinates of each Delauny triangle of the alphaShape without holes.
TriangleCoordinates_Filled = cell2mat(TriangleCoordinates_Filled); % Convert to a matrix (1 row per triangle).

HoleOrNot = ismember(TriangleCoordinates_Filled,TriangleCoordinates,'rows'); % Check row-wise which triangles are part of the point cloud (1), and which ones are part of the holes (0).
HoleTriangles = TriangleCoordinates_Filled(~HoleOrNot,:); % Extract the triangles that are part of the holes (look at the inverse of the 'HoleOrNot' vector).

% Only do this if there are holes in the alphaShape.
if ~isempty(HoleTriangles)

    % Convert the hole triangles to polyshapes, and combine the ones that
    % belong together to make up a full hole. Then do some cleaning up to
    % simplify the excess vertices/borderpoints.
    HoleTriangles = mat2cell(HoleTriangles,ones(size(HoleTriangles,1),1),size(HoleTriangles,2)); % Convert the coordinates of the hole triangles to a cell again (easier to work with).
    HoleTriangles = cellfun(@(x) reshape(x,2,3)',HoleTriangles,'UniformOutput',false); % Reshape the vector in each cell to (x,y) coordinates of the triangles.

    TrianglePolygons = cellfun(@(x) polyshape(x),HoleTriangles,'UniformOutput',false); % Convert each triangle to a polyshape.
    HolePolygons = union(vertcat(TrianglePolygons{:}),'KeepCollinearPoints',true); % Unify each triangle polyshape to make up the holes.
    HolePolygons = simplify(HolePolygons); % Simplify the unionized polyshape. This removes double vertices/borderpoints, but can potentially include holes inside the polyshape.

    IndividualHoles = num2cell(regions(HolePolygons)); % Extract a polyshape for each of the holes (can contain double coordinates or coordinates inside the polygon).
    HoleVertices = cellfun(@(x) x.Vertices,IndividualHoles,'UniformOutput',false); % Extract the vertices of these polyshapes.
    IndividualHoles = cellfun(@(x) polyshape(x),HoleVertices,'UniformOutput',false); % Make new polyshapes (this will select only the outer coordinates as vertices).

    % Extract the boundaries of each of the holes, and combine it with the
    % polyshape of the point cloud.
    [HoleCoordinates_x,HoleCoordinates_y] = cellfun(@(x) boundary(x),IndividualHoles,'UniformOutput',false); % Extract the coordinates of the borders of the holes.
    HoleCoordinates = cellfun(@(x,y) horzcat(x,y),HoleCoordinates_x,HoleCoordinates_y,'UniformOutput',false);
    BorderCoords = vertcat({BorderPoints},HoleCoordinates);

    % Make the polyshape of the point cloud.
    AllHoles = union(vertcat(IndividualHoles{:})); % Concatenate all hole polyshapes and make a single polyshape of it.
    PolyshapeCloud = subtract(PolygonCluster,AllHoles); % Remove the holes from the filled polyshape of the cluster.

% If there are no holes, then everything is just the outer boundary and
% filled polyshape.
else

    BorderCoords = {BorderPoints}; % Only the coordinates of the outer border.
    PolyshapeCloud = PolygonCluster; % Only the polyshape of the outer border shape.

end

% Turn warnings back on
warning('on','all');

end