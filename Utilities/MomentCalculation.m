function Moments = MomentCalculation(Data)
% -------------------------------------------------------------------------
% Function that calculates the scale and rotation invariant moments of a 
% data point cloud.
% Examples on how to use it:
%   PointCloudMoments = MomentCalculation(PointCloud);
% -------------------------------------------------------------------------
% Input:
%   Data:   The data cloud of which the moments should be calculated.
%           Format is a matrix containing just the (x,y) coordinates in
%           column 1 (x) and column 2 (y).
%
% Output:
%   Moments:   A structure containing all the scale invariant and rotation
%              invariant central moments
% -------------------------------------------------------------------------
% Code written by:
%   Siewert Hugelier    Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   siewert.hugelier@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   xxx
% -------------------------------------------------------------------------

% Check the input size of the data matrix
if size(Data,2) == 1
    error("The data matrix input should be specified as a 2-column matrix or as an image of size (n x m) with m > 2.")
end

if size(Data,2) == 2
    mu00 = size(Data,1);
    x = Data(:,1);
    y = Data(:,2);
    Intensity = ones(size(Data,1),1);
else
    mu00 = sum(Data(:));
    [x,y,Intensity] = find(Data);
end

Centroid(1) = sum(x.*Intensity) / mu00;
Centroid(2) = sum(y.*Intensity) / mu00;

x = x - Centroid(1);
y = y - Centroid(2);

% Calculate the scale invariant central moments of the image
Moments.n11 = sum(x.^1 .* y.^1 .* Intensity) / mu00^(1+(1+1)/2);
Moments.n20 = sum(x.^2 .* y.^0 .* Intensity) / mu00^(1+(2+0)/2);
Moments.n02 = sum(x.^0 .* y.^2 .* Intensity) / mu00^(1+(0+2)/2);
Moments.n21 = sum(x.^2 .* y.^1 .* Intensity) / mu00^(1+(2+1)/2);
Moments.n12 = sum(x.^1 .* y.^2 .* Intensity) / mu00^(1+(1+2)/2);
Moments.n30 = sum(x.^3 .* y.^0 .* Intensity) / mu00^(1+(3+0)/2);
Moments.n03 = sum(x.^0 .* y.^3 .* Intensity) / mu00^(1+(0+3)/2);

% Calculate the rotation invariant central moments of the data cloud
% This is based on (binary) image moments (e.g.: https://en.wikipedia.org/wiki/Image_moment)
Moments.M1 = Moments.n20 + Moments.n02;
Moments.M2 = (Moments.n20 - Moments.n02)^2 + 4*Moments.n11^2;
Moments.M3 = (Moments.n30 - 3*Moments.n12)^2 + (3*Moments.n21 - Moments.n03)^2;
Moments.M4 = (Moments.n30 + Moments.n12)^2 + (Moments.n21 + Moments.n03)^2;
Moments.M5 = (Moments.n30 - 3*Moments.n12)*(Moments.n30 + Moments.n12)*((Moments.n30 + Moments.n12)^2-3*(Moments.n21 + Moments.n03)^2)+(3*Moments.n21 - Moments.n03)*(Moments.n21 + Moments.n03)*(3*(Moments.n30 + Moments.n12)^2 - (Moments.n21 + Moments.n03)^2);
Moments.M6 = (Moments.n20 - Moments.n02)*((Moments.n30 + Moments.n12)^2 - (Moments.n21 + Moments.n03)^2) + 4*Moments.n11*(Moments.n30 + Moments.n12)*(Moments.n21 + Moments.n03);
Moments.M7 = (3*Moments.n21 - Moments.n03)*(Moments.n30 + Moments.n12)*((Moments.n30 + Moments.n12)^2 - 3*(Moments.n21 + Moments.n03)^2) - (Moments.n30 - 3*Moments.n12)*(Moments.n21 + Moments.n03)*(3*(Moments.n30 + Moments.n12)^2 - (Moments.n21 + Moments.n03)^2);

end