function DescriptorsTable = CalcDescriptors3D(Data)
% -------------------------------------------------------------------------
% Function that calculates the (spatial) descriptors of point clouds. Most
% descriptors/features that are calculated are calculated with respect to
% the point cloud or an accurate representation of a pointcloud (based on
% the alphaShape or polyshape functions of Matlab). A few select ones are
% calculated using a pixelized version of the point cloud.
% Example on how to use it:
%   Features = CalcDescriptors3D(PointClouds);
% -------------------------------------------------------------------------
% Input:
%   Data:   The clusters on which the descriptors should be calculated.
%           Format is a cell variable containing each cluster separately, 
%           with each cell containing just the (x,y,z) coordinates in 
%           column 1 (x), column 2 (y), and colum 3 (z). The coordinates
%           should be in nm.
%
% Output:
%   DescriptorsTable:   A table containing all the calculated descriptors.
%                       Its size will be N x 69
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

% Check if the input data is in the correct format.
if ~iscell(Data)
    error("The input should be specified as a cell array (in which each cell contains localization coordinates of a single cluster).")
end

% Force move the PLS Toolbox to the end of the path to avoid issues with
% the pca function
oldPath = path;
seperatePaths = split(oldPath,pathsep);
plsOrNot = contains(seperatePaths,'PLS_Toolbox');
newPath = join(join(seperatePaths(~plsOrNot),pathsep), join(seperatePaths(plsOrNot),pathsep));
addpath(newPath{1});

% Remove any extra columns to save memory if needed.
if size(Data{1},2) > 3
    Data = cellfun(@(x) x(:,1:3),Data,'UniformOutput',false);
end

% Determine the order of the descriptors to be calculated.
% These are divided into several groups, depending on their method of
% calculation.
GeometricDescriptors = ["NumberLocs","MajorAxis","MiddleAxis","MinorAxis","Volume","FilledVolume","ConvexVolume","GyrationRadius","EigenvaluesRatio_xy","EigenvaluesRatio_xz","EigenvaluesRatio_yz","Density","VolumeRatio","ConvexVolumeRatio","AspectRatio_xy","AspectRatio_xz","AspectRatio_yz","IsoperimetricRatio","Cubicity","Sphericity","ConvexSphericity","Convexity","EquivalentDiameter","Eigenentropy","Omnivariance","Anisotropy","Planarity","Linearity","MajorAxisNormVolume","MiddleAxisNormVolume","MinorAxisNormVolume"];
BoundaryDescriptors = ["SurfaceArea","OuterSurfaceArea","ConvexSurfaceArea","MeanCurvature"];
EllipsoidDescriptors = ["EllipsoidMajorAxis","EllipsoidMiddleAxis","EllipsoidMinorAxis","EllipsoidAspectRatio_xy","EllipsoidAspectRatio_xz","EllipsoidAspectRatio_yz","EllipsoidRatio","EquatorialEccentricity","MeridionalEccentricity"];
SkeletonDescriptors = ["SkeletonCloudWidth","SkeletonTotalLength","SkeletonIntersections","SkeletonMeanLength","SkeletonMeanOrientation_xy","SkeletonMeanOrientation_xz","SkeletonMeanOrientation_yz","SkeletonMeanTortuosity","SkeletonMeanCurvature","SkeletonTotalLengthNormVolume","SkeletonIntersectionsNormVolume"];
ImageDescriptors = ["EulerNumber","Breadth","Entropy","MeanIntensity","MeanNonZeroIntensity","RMSRoughness","Skewness"];
MomentDescriptors = "CentralMoment"+arrayfun(@string, 1:3);
FractalDescriptors = ["LocalMinkowskiBouligandDim","MinkowskiSausage","HausdorffDim","SkeletonMinkowskiSausage"];

% Combine all descriptors.
DescriptorList = [GeometricDescriptors BoundaryDescriptors EllipsoidDescriptors SkeletonDescriptors ImageDescriptors MomentDescriptors FractalDescriptors]; % Put all shape descriptor names together.

% Prepare the output table.
N = size(Data,1); % Extract the number of clusters.
DescriptorsTable = array2table(zeros(N, size(DescriptorList,2))); % Make a table with the size of the number of clusters.
DescriptorsTable.Properties.VariableNames = DescriptorList; % Give the columns of the table a proper name.

% Create a waitbar
wb = waitbar(0,['Calculating Descriptors. Please be patient... At cluster 0/' num2str(N) '     '],'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(wb,'canceling',0);

% Calculate the descriptors for each individual cluster.
for i = 1:N

    % Decide what to do when cancel is pressed.
    if getappdata(wb,'canceling')
        break
    end

    % Update the waitbar.
    waitbar(i/N,wb,['Calculating Descriptors. Please be patient... At cluster ' num2str(i) '/' num2str(N)]);

    % Extract the current cluster (easier to work with).
    Cluster = Data{i};

    % ---------------------------------------------------------------------
    % Calculate the descriptors based on the (x,y) coordinates:
    % Localization descriptors - Boundary descriptors - Dependent point
    % cloud descriptors - Ellipsoid descriptors - Humoment descriptors
    % ---------------------------------------------------------------------
    
    % Perform a Principal Component Analysis to rotate the data (maximizing
    % the length).
    [~,RotatedCluster,Eigenvalues] = pca(Cluster); % Uses the built-in Matlab pca routine.
    RotatedCluster = unique(RotatedCluster,'rows'); % This should never happen, unless really 'unlucky' with the localization. Just to avoid warning messages with alphaShape.
    Eigenvalues = Eigenvalues / sum(Eigenvalues); % Normalize the Eigenvalues between 0-1.
    x_length = max(RotatedCluster(:,1))-min(RotatedCluster(:,1)); % Determine the length of the cloud in x.
	y_length = max(RotatedCluster(:,2))-min(RotatedCluster(:,2)); % Determine the length of the cloud in y.
	z_length = max(RotatedCluster(:,3))-min(RotatedCluster(:,3)); % Determine the length of the cloud in z.
    sortedPrincipalAxes = sort([x_length y_length z_length]); % Sort the sizes from high to low.

    % Calculate the alphaShape.
    % Extract the full boundary, the outer boundary, and the boundary of 
    % the convex hull of the point cloud.
    % Extract the surface area as well.
    ClusterAlphaShape = alphaShape(RotatedCluster); % Create the alphaShape (uses Delaunay triangulation).
    AlphaOneRegion = criticalAlpha(ClusterAlphaShape,'one-region'); % Calculate the alpha corresponding to a 'one region' alphaShape, to make sure that all (x,y) coordinates are connected.
    ClusterAlphaShape.Alpha = AlphaOneRegion; % Apply it to the alphaShape.

    [I,FilledVolume] = boundary(RotatedCluster,1); % Calculate the outer border very precisely.
    OuterBorderPoints = RotatedCluster(I,:); % Extract the border coordinates.
    OuterBorderPoints = unique(OuterBorderPoints,'rows'); % Simplify of a border coordinate was taken multiple times (to avoid warning messages from alphaShape).
    ClusterOuterAlphaShape = alphaShape(OuterBorderPoints); % Create the alphaShape (uses Delaunay triangulation).
    AlphaOneRegion = criticalAlpha(ClusterOuterAlphaShape,'one-region'); % Calculate the alpha corresponding to a 'one region' alphaShape, to make sure that all (x,y) coordinates are connected.
    ClusterOuterAlphaShape.Alpha = AlphaOneRegion; % Apply it to the alphaShape.

    [I,ConvexVolume] = boundary(RotatedCluster,0.02); % Calculate the convex outer border very precisely. Put it to a small value rather than 0 to avoid too much of a simplication. Changing it from 0 to 0.02 worked really well in practise.
    ConvexBorderPoints = RotatedCluster(I,:); % Extract the border coordinates.
    ConvexBorderPoints = unique(ConvexBorderPoints,'rows'); % Simplify of a border coordinate was taken multiple times (to avoid warning messages from alphaShape).
    ClusterConvexAlphaShape = alphaShape(ConvexBorderPoints); % Create the alphaShape (uses Delaunay triangulation).
    AlphaOneRegion = criticalAlpha(ClusterConvexAlphaShape,'one-region'); % Calculate the alpha corresponding to a 'one region' alphaShape, to make sure that all (x,y) coordinates are connected.
    ClusterConvexAlphaShape.Alpha = AlphaOneRegion; % Apply it to the alphaShape.

    BorderCurvature = Curvature(OuterBorderPoints); % Only calculate the curvature for the outer boundary of the point cloud.

    % Store the shape descriptors: localization descriptors.
    DescriptorsTable(i,:).("NumberLocs") = size(Cluster,1);
    DescriptorsTable(i,:).("MajorAxis") = sortedPrincipalAxes(3);
    DescriptorsTable(i,:).("MiddleAxis") = sortedPrincipalAxes(2);
    DescriptorsTable(i,:).("MinorAxis") = sortedPrincipalAxes(1);
    DescriptorsTable(i,:).("Volume") = volume(ClusterAlphaShape);
    DescriptorsTable(i,:).("FilledVolume") = FilledVolume;
    DescriptorsTable(i,:).("ConvexVolume") = ConvexVolume;
    DescriptorsTable(i,:).("GyrationRadius") = sqrt(mean((RotatedCluster(:,1)-mean(RotatedCluster(:,1))).^2 + (RotatedCluster(:,2)-mean(RotatedCluster(:,2))).^2 + (RotatedCluster(:,3)-mean(RotatedCluster(:,3))).^2));
    DescriptorsTable(i,:).("EigenvaluesRatio_xy") = Eigenvalues(1) / Eigenvalues(2);
    DescriptorsTable(i,:).("EigenvaluesRatio_xz") = Eigenvalues(1) / Eigenvalues(3);
    DescriptorsTable(i,:).("EigenvaluesRatio_yz") = Eigenvalues(2) / Eigenvalues(3);    

    % Store the shape descriptors: boundary descriptors.
    DescriptorsTable(i,:).("SurfaceArea") = surfaceArea(ClusterAlphaShape);
    DescriptorsTable(i,:).("OuterSurfaceArea") = surfaceArea(ClusterOuterAlphaShape);
    DescriptorsTable(i,:).("ConvexSurfaceArea") = surfaceArea(ClusterConvexAlphaShape);
    DescriptorsTable(i,:).("MeanCurvature") = mean(BorderCurvature,'omitnan');

    % Store the dependent shape descriptors (coordinate-based).
    DescriptorsTable(i,:).("Density") = DescriptorsTable(i,:).("NumberLocs") / DescriptorsTable(i,:).("Volume");
    DescriptorsTable(i,:).("VolumeRatio") = DescriptorsTable(i,:).("Volume") / DescriptorsTable(i,:).("FilledVolume");
    DescriptorsTable(i,:).("ConvexVolumeRatio") = DescriptorsTable(i,:).("Volume") / DescriptorsTable(i,:).("ConvexVolume");
    DescriptorsTable(i,:).("AspectRatio_xy") = DescriptorsTable(i,:).("MajorAxis") / DescriptorsTable(i,:).("MiddleAxis");
    DescriptorsTable(i,:).("AspectRatio_xz") = DescriptorsTable(i,:).("MajorAxis") / DescriptorsTable(i,:).("MinorAxis");
    DescriptorsTable(i,:).("AspectRatio_yz") = DescriptorsTable(i,:).("MiddleAxis") / DescriptorsTable(i,:).("MinorAxis");
    DescriptorsTable(i,:).("IsoperimetricRatio") = DescriptorsTable(i,:).("SurfaceArea")^3 / DescriptorsTable(i,:).("Volume")^2;
    DescriptorsTable(i,:).("Cubicity") = DescriptorsTable(i,:).("Volume") / (DescriptorsTable(i,:).("MajorAxis")*DescriptorsTable(i,:).("MiddleAxis")*DescriptorsTable(i,:).("MinorAxis"));
    DescriptorsTable(i,:).("Sphericity") = pi^(1/3) * (6*DescriptorsTable(i,:).("FilledVolume"))^(2/3) / DescriptorsTable(i,:).("OuterSurfaceArea");
    DescriptorsTable(i,:).("ConvexSphericity") = pi^(1/3) * (6*DescriptorsTable(i,:).("ConvexVolume"))^(2/3) / DescriptorsTable(i,:).("ConvexSurfaceArea");
    DescriptorsTable(i,:).("Convexity") = DescriptorsTable(i,:).("OuterSurfaceArea") / DescriptorsTable(i,:).("ConvexSurfaceArea");
    DescriptorsTable(i,:).("EquivalentDiameter") = 2*(3/4 * DescriptorsTable(i,:).("Volume") / pi)^(1/3);
    DescriptorsTable(i,:).("Eigenentropy") = -Eigenvalues(1)*log(Eigenvalues(1)) - Eigenvalues(2)*log(Eigenvalues(2)) - Eigenvalues(3)*log(Eigenvalues(3));
    DescriptorsTable(i,:).("Omnivariance") = prod(Eigenvalues)^(1/3);
    DescriptorsTable(i,:).("Anisotropy") = (Eigenvalues(1)-Eigenvalues(3))/Eigenvalues(1);
    DescriptorsTable(i,:).("Planarity") = (Eigenvalues(2)-Eigenvalues(3))/Eigenvalues(1);
    DescriptorsTable(i,:).("Linearity") = (Eigenvalues(1)-Eigenvalues(2))/Eigenvalues(1);
    DescriptorsTable(i,:).("MajorAxisNormVolume") = DescriptorsTable(i,:).("MajorAxis") / DescriptorsTable(i,:).("Volume");
    DescriptorsTable(i,:).("MiddleAxisNormVolume") = DescriptorsTable(i,:).("MiddleAxis") / DescriptorsTable(i,:).("Volume");
    DescriptorsTable(i,:).("MinorAxisNormVolume") = DescriptorsTable(i,:).("MinorAxis") / DescriptorsTable(i,:).("Volume");

    % Calculate central moments: translation and rotational invariants. 
    % These are not scale independent, so to avoid inconsistencies
    % with image transformations, the raw coordinates (in nm) are used.
    [Moments,EllipsoidProperties] = MomentCalculation3D(RotatedCluster,CalcEllipsoid=true); % Calculates all moments in a separate function. Technically this can also be done using the rotated cluster coordinates (it should theoretically not change)! It was chosen to implement it this way to avoid numerical issues

    % Store the shape descriptors: Ellipse descriptors.
    DescriptorsTable(i,:).("EllipsoidMajorAxis") = EllipsoidProperties.MajorAxis;
    DescriptorsTable(i,:).("EllipsoidMiddleAxis") = EllipsoidProperties.MiddleAxis;
    DescriptorsTable(i,:).("EllipsoidMinorAxis") = EllipsoidProperties.MinorAxis;
    DescriptorsTable(i,:).("EllipsoidAspectRatio_xy") = DescriptorsTable(i,:).("EllipsoidMajorAxis") / DescriptorsTable(i,:).("EllipsoidMiddleAxis");
    DescriptorsTable(i,:).("EllipsoidAspectRatio_xz") = DescriptorsTable(i,:).("EllipsoidMajorAxis") / DescriptorsTable(i,:).("EllipsoidMinorAxis");
    DescriptorsTable(i,:).("EllipsoidAspectRatio_yz") = DescriptorsTable(i,:).("EllipsoidMiddleAxis") / DescriptorsTable(i,:).("EllipsoidMinorAxis");
    DescriptorsTable(i,:).("EllipsoidRatio") = 3/4 * DescriptorsTable(i,:).("Volume") / (pi * DescriptorsTable(i,:).("EllipsoidMajorAxis") * DescriptorsTable(i,:).("EllipsoidMiddleAxis") * DescriptorsTable(i,:).("EllipsoidMinorAxis") / 6);
    DescriptorsTable(i,:).("EquatorialEccentricity") = EllipsoidProperties.EquatorialEccentricity;
    DescriptorsTable(i,:).("MeridionalEccentricity") = EllipsoidProperties.MeridionalEccentricity;

    % Store the shape descriptors: 3 translation and rotation invariant 
    % Moments.
    for j = 1:3
        DescriptorsTable(i,:).("CentralMoment"+string(j)) = Moments.("J"+string(j));
    end
    
    % ---------------------------------------------------------------------
    % Calculate the descriptors based on a coordinate mask.
    % Skeleton descriptors - Image descriptors - Fractal descriptors -
    % Dependent Skeleton descriptors
    % ---------------------------------------------------------------------

    % Caclulate the image transforms needed for the next part.
    Images = ImageCloud3D(RotatedCluster); % Calculate the different images needed to calculate their proprties.
    [FractProperties, SkelProperties, GrayScaleProperties, OtherImageProperties] = ImageProperties3D(Images); % Calculate the different properties from the previously obtained images.

    % Store the shape descriptors: skeleton descriptors.
    DescriptorsTable(i,:).("SkeletonCloudWidth") = SkelProperties.Cloudwidth;
    DescriptorsTable(i,:).("SkeletonTotalLength") = SkelProperties.TotalLength;
    DescriptorsTable(i,:).("SkeletonIntersections") = SkelProperties.Intersections;
    DescriptorsTable(i,:).("SkeletonMeanLength") = SkelProperties.MeanLength;
    DescriptorsTable(i,:).("SkeletonMeanOrientation_xy") = SkelProperties.MeanOrientation_xy;
    DescriptorsTable(i,:).("SkeletonMeanOrientation_xz") = SkelProperties.MeanOrientation_xz;
    DescriptorsTable(i,:).("SkeletonMeanOrientation_yz") = SkelProperties.MeanOrientation_yz;
    DescriptorsTable(i,:).("SkeletonMeanTortuosity") = SkelProperties.MeanTortuosity;
    DescriptorsTable(i,:).("SkeletonMeanCurvature") = SkelProperties.MeanCurvature;    

    % Store the shape descriptors: Image descriptors.
    DescriptorsTable(i,:).("Entropy") = GrayScaleProperties.Entropy;
    DescriptorsTable(i,:).("MeanIntensity") = GrayScaleProperties.MeanIntensity;
    DescriptorsTable(i,:).("MeanNonZeroIntensity") = GrayScaleProperties.MeanNonZeroIntensity;
    DescriptorsTable(i,:).("RMSRoughness") = GrayScaleProperties.RMSRoughness;
    DescriptorsTable(i,:).("Skewness") = GrayScaleProperties.Skewness;
    DescriptorsTable(i,:).("EulerNumber") = OtherImageProperties.EulerNumber;
    DescriptorsTable(i,:).("Breadth") = OtherImageProperties.Breadth;

    % Store the shape descriptors: fractal descriptors.
    DescriptorsTable(i,:).("LocalMinkowskiBouligandDim") = FractProperties.MBDim;
    DescriptorsTable(i,:).("MinkowskiSausage") = FractProperties.MSausageDim;
    DescriptorsTable(i,:).("HausdorffDim") = FractProperties.HDim;
    DescriptorsTable(i,:).("SkeletonMinkowskiSausage") = FractProperties.MSausageDim_Skeleton;

    % Store the dependent shape descriptors.
    DescriptorsTable(i,:).("SkeletonTotalLengthNormVolume") = DescriptorsTable(i,:).("SkeletonTotalLength") ./  DescriptorsTable(i,:).("Volume");
    DescriptorsTable(i,:).("SkeletonIntersectionsNormVolume") = DescriptorsTable(i,:).("SkeletonIntersections") ./  DescriptorsTable(i,:).("Volume");

end

% Update the waitbar and then close it.
waitbar(1,wb,'Calculating Descriptors. Finished...');
delete(wb);

% Turn warnings back on
warning('on','all');

end