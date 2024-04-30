function Images = ImageCloud3D(Pointcloud)
% -------------------------------------------------------------------------
% Function that is the umbrella function to calculate different image
% transforms. It makes different image representations of the coordinates 
% of the point cloud. To be consistent between structures, it is advised to
% perform a PCA rotation of the coordinates beforehand to set the x-axis as 
% the major variance direction.
% Examples on how to use it:
%   [~,RotatedPointcloud] = pca(Pointcloud);
%   Images = ImageCloud3D(RotatedPointcloud);
% -------------------------------------------------------------------------
% Input:
%   Pointcloud:     The point cloud data (this is a n x 3 matrix).
%
% Output:
%   Images:             The image transforms of the point cloud data (i.e.
%                       a structure containing the logical, boundary, 
%                       Euclidean distance, skeleton and branchpoint
%                       images).
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
if size(Pointcloud,2) ~= 3
    error("The input should be specified as a n x 3 matrix.");
end

% Resample the cluster to match the axial precision. In our data, this is
% ~50nm.
Scale = 20; % In nanometers.
Pointcloud = (Pointcloud - min(Pointcloud))/Scale + 1; % We add the 1 to avoid issues while making the image.

% Extract the image sizes. It should not be larger than the biggest values
% in our point cloud. Also make an empty 3D image to construct our images.
ImgSize = ceil(max(Pointcloud));
GrayScaleImage = zeros(ImgSize(2),ImgSize(1),ImgSize(3)); % x and y are inversed because rows should be y, columns should be x.

% Make the GrayScale image of the coordinates. A voxel here is considered
% to have the size of the axial precision (1 x 1 x 1 pixel = 50 x 50 x 50 
% nm).
for i = 1:size(Pointcloud,1)
    GrayScaleImage(round(Pointcloud(i,2)),round(Pointcloud(i,1)),round(Pointcloud(i,3))) = GrayScaleImage(round(Pointcloud(i,2)),round(Pointcloud(i,1)),round(Pointcloud(i,3))) + 1;
end

% Convert the grayscale image to a logical image and clean it up to
% facilitate further image transformations.
LogicalImage = logical(GrayScaleImage);
LogicalImage = imclose(LogicalImage,strel("sphere",3)); % Close the image to avoid having disconnected pixels.
LogicalImage = bwmorph3(LogicalImage,'clean'); % Clean the image from isolated voxels.

% Transform the logical coordinate mask into different other images.
EuclideanDistImage = bwdist(~LogicalImage); % Calculate the Euclidean distance image (closest distance to a border pixel).
SkeletonImage = bwskel(LogicalImage);
BranchPointImage = bwmorph3(SkeletonImage, 'branchpoints');
BoundaryImage = bwperim(LogicalImage,6);

% Add the different images to the output variable (a structure).
Images.GrayScaleImage = GrayScaleImage; % Grayscale Image of the pointcloud.
Images.Logical = LogicalImage; % Logical image of the pointcloud.
Images.BoundaryImage = BoundaryImage; % Boundary image of the pointcloud.
Images.EuclideanDistImage = EuclideanDistImage; % Euclidean distance image of the pointcloud.
Images.SkeletonImage = SkeletonImage; % Skeleton image of the pointcloud.
Images.BranchPointImage = BranchPointImage; % Branchpoint image of the pointcloud.

% The conversion of pixels in the images to size in the original point
% cloud. 
Images.Scale = Scale; % The scale of the voxels related to the point cloud data in nm.

end