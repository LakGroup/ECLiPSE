function ImageTransforms = ImageTransformCloud(PolygonCloud)
% -------------------------------------------------------------------------
% Function that is the umbrella function to calculate different image
% transforms. It makes an accurate representation of a polyshape of the
% point cloud and then calculates the different images from that.
% Examples on how to use it:
%   PointCloud_PS = polyshape(PointCloud);
%   Images = ImageTransformCloud(PointCloud_PS);
% -------------------------------------------------------------------------
% Input:
%   PolygonCloud:   The polygon of the point cloud data (this is a
%                   polyshape).
%
% Output:
%   ImageTransforms:    The image transforms of the point cloud data (i.e.
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
if ~isa(PolygonCloud,'polyshape')
    error("The input should be specified as a polyshape (see 'help polyshape')");
end

% Extract the length scales at which the image should be displayed.
minx = min(PolygonCloud.Vertices(:,1)); % Extract the minimum in the x-direction.
maxx = max(PolygonCloud.Vertices(:,1)); % Extract the maximum in the x-direction.
miny = min(PolygonCloud.Vertices(:,2)); % Extract the minimum in the y-direction.
maxy = max(PolygonCloud.Vertices(:,2)); % Extract the maximum in the y-direction.

% Plot the polyshape of the point cloud and convert it into a highly 
% detailed image (the point cloud will be upscales depending on the size of
% it).
% To do this, we use the 'capture frame' option of MATLAB and then convert
% that frame into a logical image.
fig = figure('Visible','off','MenuBar','none'); % Make a new figure but make it invisible.
set(fig,'color','w'); % Set the background of the figure to white (for easier control).
plot(PolygonCloud,'FaceColor','k','FaceAlpha',1); % Plot the polyshape of the point cloud in black (no transparency).
axis([minx maxx miny maxy]); % Set the axes so that it fills the entire figure.
axis off; % Turn off the axis (it does not have to be printed).
axis equal;

CapturedFrame = getframe(fig); % Extract the plotted frame.
close(fig); % close the figure.
RGB_Image = frame2im(CapturedFrame); % Extract the data from the plotted frame (which is now an RGB image).
Red_Image = RGB_Image(:,:,1); % Only extract the first colour (the rest is just repetition for our sake).
Red_Image(Red_Image~=255) = 1; % Set all values that are not 'white' (w = 255 here) to 1 (i.e., part of the point cloud).
Red_Image(Red_Image==255) = 0; % Set the remaining values to 0 (i.e., not part of the point cloud).

% Clean up the image as there is now a lot of blank space (as the entire
% figure got captured, and not just the point cloud. The image is then
% cropped to those dimensions, and afterwards padded by 3 pixels.
[IdxRow,IdxCol] = find(Red_Image); % Find the indices of the point cloud.
Red_Image = Red_Image(min(IdxRow):max(IdxRow),min(IdxCol):max(IdxCol)); % Crop the image to only the point cloud.
CoordinateImage = logical(padarray(double(Red_Image), [3 3], 0, "both")); % Pad the image by 3 pixels (arbitrary choice) to avoid problems with the skeletonization.

% Transform the logical coordinate mask into different other images.
BoundaryImage = bwperim(CoordinateImage); % Make the image containing all the boundary pixels.
EuclideanDistImage = bwdist(~CoordinateImage); % Calculate the Euclidean distance image (closest distance to a border pixel).
SkeletonImage = bwmorph(CoordinateImage,'skel',Inf); % Make the skeleton image of the point cloud.
SkeletonImage = bwmorph(SkeletonImage,'thin',Inf); % Thin the skeleton image because sometimes branchpoints are awkward in the skeletonized image.
% SkeletonImage = bwskel(CoordinateImage); % Make the skeleton image of the point cloud.
BranchPointImage = bwmorph(SkeletonImage, 'branchpoints',Inf); % Extract the branchpoints of the skeletonized image.

% Add the different images to the output variable (a structure).
ImageTransforms.Logical = CoordinateImage; % Logical image of the coordinates.
ImageTransforms.Boundary = BoundaryImage; % Image of the boundary pixels.
ImageTransforms.EucDist = double(EuclideanDistImage); % Euclidean distance image.
ImageTransforms.Skeleton = SkeletonImage; % Skeletonized image.
ImageTransforms.BranchPoint = BranchPointImage; % Branchpoints of the skeleton.

% The conversion of pixels in the new images to size in the original image. 
% The size may be different for x and y because of rounding errors.
% On our microscope (1 px = 117 nm), the difference between x/y values for 
% 12 000 clusters of varying size (experimental data):
%   - Min: 0 nm; Max: 0.7826 nm
%   - Median: 0.0031 nm; Mean: 0.0062 nm
% We deemed these values negligible and will take the mean of x/y in future
% calculations.
ImageTransforms.Scale = [(maxy-miny)/(size(BoundaryImage,1)-6) (maxx-minx)/(size(BoundaryImage,2)-6)]; % The scale of the point cloud data in pixels (needed to convert back to true distances, -6 because there is padding of 3 pixels in each direction). Keep in mind that x is columns (2nd MATLAB direction) and y is rows (1st MATLAB direction).

end