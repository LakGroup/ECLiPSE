function [FractProperties, SkelProperties] = ImageProperties(Image)
% -------------------------------------------------------------------------
% Function that calculates a range of image properties (fractals and
% skeleton properties). The images in this work are accurate
% representations of the point cloud (based on alphaShape / polyshape), but
% this function can be used for any type of image.
% Examples on how to use it:
%   PointCloud_PS = polyshape(PointCloud);
%   Images = ImageTransformCloud(PointCloud_PS);
%   [Fractals,Skeleton] = ImageProperties(Images);
% -------------------------------------------------------------------------
% Input:
%   Image:   A structure containing different image tranforms (coming from
%            the 'ImageTransformCloud' function).
%
% Output:
%   FractProperties:    Different calculated fractal properties of the 
%                       point cloud data using the image transforms.
% -------------------------------------------------------------------------
% Code written by:
%   Siewert Hugelier    Lakadamyali lab, University of Pennsylvania (USA)
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
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
if ~isstruct(Image)
    error("The input should be specified as a structure containing 5 different images: (i) Logical image, (ii) Boundary image, (iii) Euclidean Distance image, (iv) Skeleton image, and (v) Branchpoint image. (see 'help ImageTransformCloud' for more info)");
end

% Extract the different images as it is easier to work with.
CoordinateImage = Image.Logical; % Logical image of the coordinates.
BoundaryImage = Image.Boundary; % Image of the boundary pixels.
EuclideanDistImage = Image.EucDist; % Euclidean distance image.
SkeletonImage = Image.Skeleton; % Skeletonized image.
Scale = mean(Image.Scale); % The conversion of pixels in the new images to size in the original image. The size may be different for x and y because of rounding errors.

% Extract the individual lines/branches from the skeleton image.
[Branches, NumberBranchPoints] = ExtractBranches(SkeletonImage);

% Calculate fractal properties based on these obtained images, and add them
% to an output variable (a structure).
FractProperties.MBDim = LocalMinkowskiBouligandDim(CoordinateImage); % Calculate the Minkowski-Bouligand (box counting) dimension of the point cloud.
FractProperties.MSausageDim = MinkowskiSausageDim(CoordinateImage); % Calculate the Minkowski curve of the point cloud.
FractProperties.MSausageDim_Skeleton = MinkowskiSausageDim(SkeletonImage); % Calculate the Minkowski curve of the skeleton image of the point cloud.
FractProperties.HDim = HausdorffDim(CoordinateImage); % Calculate the Haussdorff dimension of the point cloud (i.e., a measure of the roughness).
FractProperties.HDim_Boundary = HausdorffDim(BoundaryImage); % Calculate the Haussdorff dimension of the border of the point cloud (i.e., a measure of the roughness).
FractProperties.HDim_Skeleton = HausdorffDim(SkeletonImage); % Calculate the Haussdorff dimension of the skeleton image of the point cloud (i.e., a measure of the roughness).

% Calculations based on the skeleton image, and add them to an output
% variable (a structure).
SkeletonRadius = EuclideanDistImage(SkeletonImage); % Extract the shortest distance to the border of the skeleton image.
BranchAngles = cellfun(@(x) atan2(x(1,2) - x(end,2),x(1,1) - x(end,1)),Branches); % Calculate the angle between the first and last point of the branch.
BranchDistance = cellfun(@(x) pdist2(x(1,:),x(end,:)),Branches); % Calculate the distance covered by the branch.
BranchLength = cellfun(@(x) size(x,1),Branches);
SkeletonCurvature = mean(cellfun(@(x) mean(x,'omitnan'),cellfun(@(x) Curvature(x(:,1:2)),Branches,'UniformOutput',false)),'omitnan');

SkelProperties.Cloudwidth = median(SkeletonRadius)*2 * Scale; % The median width of the point cloud (expressed in pixels of the original image).
SkelProperties.TotalLength = sum(SkeletonImage(:)) * Scale; % The total length of the skeleton image (expressed in pixels of the original image).
SkelProperties.Intersections = NumberBranchPoints; % Save the number of intersections (i.e., branchpoints).
SkelProperties.MeanLength = mean(BranchLength) * Scale; % The mean length of the branches (expressed in pixels of the original image).
SkelProperties.MeanOrientation = mean(abs(rad2deg(BranchAngles))); % The mean orientation of the branches (in degrees).
SkelProperties.MeanTortuosity = mean(BranchLength./BranchDistance); % The mean tortuosity of the branches.
SkelProperties.MeanCurvature = SkeletonCurvature; % The mean curvature of the branches.

end