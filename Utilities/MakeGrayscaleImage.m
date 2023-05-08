function [GrayscaleImage,GrayscaleProperties] = MakeGrayscaleImage(DataCloud, UpscaleFactor, PadFactor, DilationFactor)
% -------------------------------------------------------------------------
% Function that makes a grayscale image from a data point cloud. It uses an
% upscale factor to supersample the pixelized image to get more accurate
% representations of the point cloud (continuous coordinates).
% Examples on how to use it:
%   [~,GrayscaleProperties] = MakeGrayscaleImage(PointCloud, 10, 3, 1);
% -------------------------------------------------------------------------
% Input:
%   DataCloud:      The data cloud of which the grayscale image should be 
%                   calculated. Format is a matrix containing just the 
%                   (x,y) coordinates in column 1 (x) and column 2 (y).
%   UpscaleFactor:  The factor with which the coordinates should be
%                   upscales to preserve details in the gray scale image.
%   Padfactor:      The factor with which the grayscale image should be
%                   padded before perfoming image transformation.
%   DilationFactor: The factor with which the grayscale image will be
%                   eroded and dilated to preserve details in its shapes.
%
% Output:
%   GrayscaleImage:         The grayscale image of the data cloud. This 
%                           will be an image that has been upscaled to get 
%                           more detail.
%   GrayscaleProperties:    Soem texture properties from the gray
% -------------------------------------------------------------------------
% Code written by:
%   Siewert Hugelier    Lakadamyali lab, University of Pennsylvania (USA)
%   Hannah Kim          Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   siewert.hugelier@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   xxx
% -------------------------------------------------------------------------

% Check the input size of the data matrix
if size(DataCloud,2) > 2 || size(DataCloud,2) == 1
    error("The data matrix input should be specified as a 2-column matrix.")
end

if numel(UpscaleFactor) ~= 1 || numel(PadFactor) ~= 1 || numel(DilationFactor) ~= 1
    error("The upscale factor, padding factor and dilation factor should be provided as a number.")
end

% Upscale the data. This will spread out coordinates within the same pixel
% over multiple pixels and preserve some of the cloud details.
DataCloud = DataCloud * UpscaleFactor; % Multiply the data coordinated with the upscale factor.

% Transform continuous coordinates to image pixels, and make a logical mask
% of active and inactive pixels.
minx = min(DataCloud(:,1)); % Minimum value in the x-direction.
maxx = max(DataCloud(:,1)); % Maximum value in the x-direction.
miny = min(DataCloud(:,2)); % Minimum value in the y-direction.
maxy = max(DataCloud(:,2)); % Maximum value in the y-direction.
xdim = ceil(maxx - minx); % Dimension of the cloud in x (columns of the images).
ydim = ceil(maxy - miny); % Dimension of the cloud in y (rows of the images).

Pixel_idx = sub2ind([ydim xdim], floor(DataCloud(:,2)-miny+1),floor(DataCloud(:,1)-minx+1)); % Transform the coordinates to image pixels (pixels start with 1, hence the +1).
LogicalMask = false([ydim xdim]); % Make an empty logical mask of the coordinates.
LogicalMask(Pixel_idx) = 1; % Set the data cloud coordinates to 1.

LogicalMask = padarray(LogicalMask, [PadFactor PadFactor], 0, 'both'); % Pad the image to not affect the dilation/erosion.
LogicalMask = imdilate(LogicalMask,strel('disk',DilationFactor)); % Dilate the image to fill holes (& remove single pixels).
LogicalMask = imerode(LogicalMask,strel('disk',DilationFactor)); % Erode the image to remove the excess area pixels.

% Make the actual grayscale image by using the 3D histogram function and
% then smoothing out the outliers (they can be due to artefacts etc.).
GrayscaleImage = hist3(DataCloud, [xdim ydim])';
GrayscaleImage = padarray(GrayscaleImage, [PadFactor PadFactor], 0, 'both'); % Pad the image to not affect the dilation/erosion.
GrayscaleImage = imdilate(GrayscaleImage,strel('disk',DilationFactor)); % Dilate the image to fill holes (& remove single pixels).
GrayscaleImage = imerode(GrayscaleImage,strel('disk',DilationFactor)); % Erode the image to remove the excess area pixels.

UniqueIntensities = unique(GrayscaleImage(:)); % Get the unique intensity values of the grayscale image.
Outliers = UniqueIntensities(isoutlier(UniqueIntensities,'quartiles')); % Determine whether or not there are any otuliers in these values.

if ~isempty(Outliers)
    [OutlierId_Row, OutlierId_Col] = find(ismember(GrayscaleImage,Outliers)); % Find the pixel coordinates of the outliers
    OutlierIDs = sub2ind(size(GrayscaleImage),OutlierId_Row,OutlierId_Col); % Transform the coordinates of the outlier pixels to pixel IDs
    
    GrayScaleImage_Conv = conv2(GrayscaleImage,[0.5 1 0.5; 1 0 1; 0.5 1 0.5],'same') / 6; % Calculate the weighted average of the neighbouring pixels.
    GrayscaleImage(OutlierIDs) = GrayScaleImage_Conv(OutlierIDs); % Replace the grayscale image values by these averages.
end

% Remove any pixels that should not contain any values due to some odd
% image transformations, then crop the image again to its original
% dimensions.
GrayscaleImage = GrayscaleImage .* LogicalMask; % Multiply the grayscale image with the logical image.
GrayscaleImage = GrayscaleImage(PadFactor+1:end-PadFactor,PadFactor+1:end-PadFactor); % Remove the padding from the image.

% Calculate grayscale image properties, based on the co-occurence matrix
Comatrix = graycomatrix(GrayscaleImage); % Calculate the gray-level co-occurence matrix.
CoProps = graycoprops(Comatrix); % Calculate the properties of this gray-level co-occurence matrix.
GrayscaleProperties.Contrast = CoProps.Contrast; % A measure of the intensity contrast between pixels and their neighbours.
GrayscaleProperties.Correlation = CoProps.Correlation; % A measure of the correlation between neighbouring pixels.
GrayscaleProperties.Energy = CoProps.Energy; % A measure for the uniformity of the image.

% Calculate other grayscale image properties.
GrayscaleVectorized = GrayscaleImage(:); % Vectorized form of the grayscale image.

GrayscaleProperties.Entropy = entropy(GrayscaleImage); % A measure of randomness to characterize the texture of the image.
GrayscaleProperties.Kurtosis = kurtosis(GrayscaleVectorized); % The kurtosis of the grayscale image.
GrayscaleProperties.MeanIntensity = mean(GrayscaleVectorized); % Mean intensity of the grayscale image.
GrayscaleProperties.MeanNonZeroIntensity = mean(nonzeros(GrayscaleVectorized)); % Mean non-zero intensity of the grayscale image.
GrayscaleProperties.RMSRoughness = sqrt(mean(GrayscaleVectorized.^2)); % The RMS roughness of the grayscale image.
GrayscaleProperties.Skewness = skewness(GrayscaleVectorized); % The skewness of the grayscale image.

end
