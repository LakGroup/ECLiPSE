function [FractProperties, SkelProperties, GrayScaleProperties, OtherImageProperties] = ImageProperties3D(Images)
% -------------------------------------------------------------------------
% Function that calculates a range of image properties (fractals, skeleton,
% Grayscale and other image properties).
% Examples on how to use it:
%   Images = ImageCloud3D(Pointcloud)
%   [Fractals, Skeleton, GrayScale, Other] = ImageProperties3D(Images);
% -------------------------------------------------------------------------
% Input:
%   Images:   A structure containing different image tranforms (coming from
%            the 'ImageCloud3D' function).
%
% Output:
%   FractProperties:        Different calculated fractal properties of the 
%                           point cloud data using the image transforms.
%   SkelProperties:         Different calculated skeleton properties of the 
%                           point cloud data using the image transforms.
%   GrayScaleProperties:    Different calculated grayscale properties of 
%                           the point cloud data using the image 
%                           transforms.
%   OtherImageProperties:   Different calculated other image properties 
%                           of the point cloud data using the image 
%                           transforms.
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
if ~isstruct(Images)
    error("The input should be specified as a structure containing at least 4 different images: (i) GrayScale image, (ii) Logical image, (iii) Euclidean Distance image, (iv) Skeleton image. It should also contain the scale of the voxels (in nm). (see 'help ImageCloud3D' for more info)");
end

% Extract the different images as it is easier to work with.
GrayScaleImage = Images.GrayScaleImage; % Grayscale Image of the pointcloud.
LogicalImage = Images.Logical; % Logical image of the pointcloud.
EuclideanDistImage = Images.EuclideanDistImage; % Euclidean distance image of the pointcloud.
SkeletonImage = Images.SkeletonImage; % Skeleton image of the pointcloud.
Scale = Images.Scale; % The conversion of voxels in the new images to size in the original image.

% [m, n, o] = size(SkeletonImage);
% BP_Image = find(SkeletonImage(:));
% [bp_x,bp_y,bp_z] = ind2sub([m, n, o],BP_Image);
% figure(1);
% plot3(bp_x,bp_y,bp_z,'.k','MarkerSize',10);
% drawnow
% pause(0.5)

% Extract the individual lines/branches from the skeleton image.
[Branches, NumberBranchPoints] = ExtractBranches3D(SkeletonImage);
BranchLengths = cellfun('size',Branches,1)<3;
Branches(BranchLengths) = [];

% Calculate fractal properties based on these obtained images, and add them
% to an output variable (a structure).
MBDim = LocalMinkowskiBouligandDim(LogicalImage); % Calculate the Minkowski-Bouligand (box counting) dimension of the point cloud.
if isnan(MBDim)
    MBDim = 0;
end
HDim = HausdorffDim(LogicalImage);
if isnan(HDim)
    HDim = 0;
end
MSausageDim = MinkowskiSausageDim(LogicalImage);
if isnan(MSausageDim)
    MSausageDim = 0;
end
MSausageDim_Skeleton = MinkowskiSausageDim(SkeletonImage);
if isnan(MSausageDim_Skeleton)
    MSausageDim_Skeleton = 0;
end
FractProperties.MBDim = MBDim; % Store the Minkowski-Bouligand (box counting) dimension of the point cloud.
FractProperties.MSausageDim = MSausageDim; % Calculate the Minkowski curve of the point cloud.
FractProperties.MSausageDim_Skeleton = MSausageDim_Skeleton; % Calculate the Minkowski curve of the skeleton image of the point cloud.
FractProperties.HDim = HDim; % Calculate the Haussdorff dimension of the point cloud (i.e., a measure of the roughness).

% Calculations based on the skeleton image, and add them to an output
% variable (a structure).
SkeletonRadius = EuclideanDistImage(SkeletonImage); % Extract the shortest distance to the border of the skeleton image.

if ~isempty(Branches)
    BranchAngles_xy = cellfun(@(x) atan2(x(1,2) - x(end,2),x(1,1) - x(end,1)),Branches); % Calculate the angle between the first and last point of the branch in x-y.
    BranchAngles_xz = cellfun(@(x) atan2(x(1,3) - x(end,3),x(1,1) - x(end,1)),Branches); % Calculate the angle between the first and last point of the branch in x-z.
    BranchAngles_yz = cellfun(@(x) atan2(x(1,3) - x(end,3),x(1,2) - x(end,2)),Branches); % Calculate the angle between the first and last point of the branch in y-z.
    BranchDistance = cellfun(@(x) pdist2(x(1,:),x(end,:)),Branches); % Calculate the distance covered by the branch.
else
    BranchAngles_xy = 0;
    BranchAngles_xz = 0;
    BranchAngles_yz = 0;
    BranchDistance = 0;
end

BranchLength = cellfun(@(x) size(x,1),Branches); % Calculate the length covered by the branch.
if isempty(BranchLength)
    BranchLength = 0;
end

Tortuosity = BranchLength./BranchDistance;
Tortuosity = Tortuosity(~isinf(Tortuosity));

if ~isempty(Branches)
    SkeletonCurvature = mean(cellfun(@(x) mean(x,'omitnan'),cellfun(@(x) Curvature(x(:,1:3)),Branches,'UniformOutput',false)),'omitnan');
else
    SkeletonCurvature = 0;
end

SkelProperties.Cloudwidth = median(SkeletonRadius)*2 * Scale; % The median width of the point cloud (expressed in voxels of the original image).
SkelProperties.TotalLength = sum(SkeletonImage(:)) * Scale; % The total length of the skeleton image (expressed in voxels of the original image).
SkelProperties.Intersections = NumberBranchPoints; % Save the number of intersections (i.e., branchpoints).
SkelProperties.MeanLength = mean(BranchLength,'omitnan') * Scale; % The mean length of the branches (expressed in voxels of the original image).
SkelProperties.MeanOrientation_xy = mean(abs(rad2deg(BranchAngles_xy)),'omitnan'); % The mean orientation of the branches (in degrees).
SkelProperties.MeanOrientation_xz = mean(abs(rad2deg(BranchAngles_xz)),'omitnan'); % The mean orientation of the branches (in degrees).
SkelProperties.MeanOrientation_yz = mean(abs(rad2deg(BranchAngles_yz)),'omitnan'); % The mean orientation of the branches (in degrees).
SkelProperties.MeanTortuosity = mean(Tortuosity,'omitnan'); % The mean tortuosity of the branches.
if isnan(SkelProperties.MeanTortuosity)
    SkelProperties.MeanTortuosity = 0;
end
SkelProperties.MeanCurvature = SkeletonCurvature; % The mean curvature of the branches.

% Calculate grayscale image properties.
% Some properties are calculated using the matImage toolbox functions
% (https://github.com/mattools/matImage). The functions used are included 
% in this analysis package as well to avoid any issues with compatibility.
GrayscaleVectorized = GrayScaleImage(:); % Vectorized form of the grayscale image.

GrayScaleProperties.Entropy = imEntropy(GrayScaleImage);
GrayScaleProperties.MeanIntensity = mean(GrayscaleVectorized,'omitnan'); % Mean intensity of the grayscale image.
GrayScaleProperties.MeanNonZeroIntensity = mean(nonzeros(GrayscaleVectorized),'omitnan'); % Mean non-zero intensity of the grayscale image.
GrayScaleProperties.RMSRoughness = sqrt(mean(GrayscaleVectorized.^2,'omitnan')); % The RMS roughness of the grayscale image.
GrayScaleProperties.Skewness = skewness(GrayscaleVectorized); % The skewness of the grayscale image.

% Calculate the properties based of other image. These properties are based
% of the matImage toolbox (https://github.com/mattools/matImage). The
% functions used are included in this analysis package as well to avoid any
% issues with compatibility.
OtherImageProperties.EulerNumber = imEuler3d(LogicalImage);
OtherImageProperties.Breadth = imMeanBreadth(LogicalImage) * Scale;

end