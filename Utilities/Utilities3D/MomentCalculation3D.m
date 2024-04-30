function [Moments,EllipsoidProperties] = MomentCalculation3D(Data,varargin)
% -------------------------------------------------------------------------
% Function that calculates the scale and rotation invariant moments of a 
% data point cloud.
% Examples on how to use it:
%   PointCloudMoments = MomentCalculation(PointCloud);
% -------------------------------------------------------------------------
% Input:
%   Data:   The data cloud of which the moments should be calculated.
%           Format is a matrix containing just the (x,y,z) coordinates in
%           column 1 (x), column 2 (y), and colum 3 (z).
%           Alternative format can also be a m x n x o binary image.
%
% Optional Input:
%   CalcEllipsoid:  A 'true' or 'false' flag to calculate the ellipsoid
%                   parameters.
%
% Output:
%   Moments:    A structure containing all the translation invariant and 
%               rotation invariant central moments
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

% Parse the input.
DefaultCalcEllipsoid = true;
p = inputParser;
addRequired(p,'Data');
addOptional(p,'CalcEllipsoid',DefaultCalcEllipsoid,@islogical);
parse(p,Data,varargin{:});

Data = p.Results.Data;
CalcEllipsoid = p.Results.CalcEllipsoid;

% Check the input size of the data matrix.
if ismatrix(Data) && size(Data,2) ~= 3
    error("The data matrix input should be specified as a 3-column matrix or as a 3-dimensional image.")
end
if ~ismatrix(Data) && ndims(Data) ~= 3
    error("The data matrix input should be specified as a 3-column matrix or as a 3-dimensional image.")
end

% Transform the data to calculate the moments around the centroid.
if size(Data,2) == 3

    % Prepare the data.
    mu000 = size(Data,1);
    x = Data(:,1);
    y = Data(:,2);
    z = Data(:,3);

    Centroid(1) = sum(x) / mu000;
    Centroid(2) = sum(y) / mu000;
    Centroid(3) = sum(z) / mu000;

    x = x - Centroid(1);
    y = y - Centroid(2);
    z = z - Centroid(3);

    % Calculate the scale invariant central moments of the image.
    % Formulas from: https://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/FISHER/MOMINV3D/inv3d.htm
    Moments.mu200 = sum(x.^2 .* y.^0 .* z.^0) / mu000 + 1/12;
    Moments.mu110 = sum(x.^1 .* y.^1 .* z.^0) / mu000;
    Moments.mu101 = sum(x.^1 .* y.^0 .* z.^1) / mu000;
    Moments.mu020 = sum(x.^0 .* y.^2 .* z.^0) / mu000 + 1/12;
    Moments.mu011 = sum(x.^0 .* y.^1 .* z.^1) / mu000;
    Moments.mu002 = sum(x.^0 .* y.^0 .* z.^2) / mu000 + 1/12;
else

    % Prepare the data.
    [m,n,o] = size(Data);
    [r,c,v] = ind2sub([m n o],find(Data));
    Centroid = mean([c r v],1);

    x = ((1:n)-Centroid(1));
    y = ((1:m)-Centroid(2))';
    z = (1:o)-Centroid(3);

    % Calculate the scale invariant central moments of the image.
    mu000 = sum(Data(:));
    mu200 = (y.^2 * x.^0) .* reshape(z.^0,[1 1 o]) .* Data;
    Moments.mu200 = sum(mu200(:)) / mu000 + 1/12;
    mu110 = (y.^1 * x.^1) .* reshape(z.^0,[1 1 o]) .* Data;
    Moments.mu110 = sum(mu110(:)) / mu000;
    mu101 = (y.^1 * x.^0) .* reshape(z.^1,[1 1 o]) .* Data;
    Moments.mu101 = sum(mu101(:)) / mu000;
    mu020 = (y.^0 * x.^2) .* reshape(z.^0,[1 1 o]) .* Data;
    Moments.mu020 = sum(mu020(:)) / mu000 + 1/12;
    mu011 = (y.^0 * x.^1) .* reshape(z.^1,[1 1 o]) .* Data;
    Moments.mu011 = sum(mu011(:)) / mu000;
    mu002 = (y.^0 * x.^0) .* reshape(z.^2,[1 1 o]) .* Data;
    Moments.mu002 = sum(mu002(:)) / mu000 + 1/12;

end

% Calculate the translation/rotation invariant central moments of the data
% cloud.
% Formulas from: https://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/FISHER/MOMINV3D/inv3d.htm
Moments.J1 = Moments.mu200 + Moments.mu020 + Moments.mu002;
Moments.J2 = Moments.mu200*Moments.mu020 + Moments.mu200*Moments.mu002 + Moments.mu020*Moments.mu002 - Moments.mu110^2 - Moments.mu101^2 - Moments.mu011^2;
Moments.J3 = Moments.mu200*Moments.mu020*Moments.mu002 + 2*Moments.mu110*Moments.mu101*Moments.mu011 - Moments.mu002*Moments.mu110^2 - Moments.mu020*Moments.mu101^2 - Moments.mu200*Moments.mu011^2;

% Calculate the properties of the ellipsoid that has the same second
% central moment as the data.
if CalcEllipsoid == 1

    % Construct the covariance matrix of the central moments and calculate
    % the eigenvalues (sorted from high to low).
    CovMatrix = [Moments.mu200 Moments.mu110 Moments.mu101;Moments.mu110 Moments.mu020 Moments.mu011;Moments.mu101 Moments.mu011 Moments.mu002];
    [~,Eigenvalues] = svd(CovMatrix);
    Eigenvalues = sort(diag(Eigenvalues), 'descend');

    % Calculate the ellipsoid axes
    EllipsoidProperties.MajorAxis = 2 * sqrt(5 * Eigenvalues(1));
    EllipsoidProperties.MiddleAxis = 2 * sqrt(5 * Eigenvalues(2));
    EllipsoidProperties.MinorAxis = 2 * sqrt(5 * Eigenvalues(3));

    % Calculate the eccentricity of the ellipsoid. The formula is 
    % simplified from the 2D code.
    EllipsoidProperties.EquatorialEccentricity = sqrt(1 - (Eigenvalues(3) / Eigenvalues(2))^2);
    EllipsoidProperties.MeridionalEccentricity = sqrt(1 - (Eigenvalues(3) / Eigenvalues(1))^2);

else

    EllipsoidProperties = [];

end