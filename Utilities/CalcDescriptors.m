function DescriptorsTable = CalcDescriptors(Data)
% -------------------------------------------------------------------------
% Function that calculates the (spatial) descriptors of point clouds. Most
% descriptors/features that are calculated are calculated with respect to
% the point cloud or an accurate representation of a pointcloud (based on
% the alphaShape or polyshape functions of Matlab). A few select ones are
% calculated using a pixelized version of the point cloud.
% Example on how to use it:
%   Features = CalcDescriptors(PointClouds);
% -------------------------------------------------------------------------
% Input:
%   Data:   The clusters on which the descriptors should be calculated.
%           Format is a cell variable containing each cluster separately, 
%           with each cell containing just the (x,y) coordinates in column 
%           1 (x) and column 2 (y).
%
% Output:
%   DescriptorsTable:   A table containing all the calculated descriptors.
%                       Its size will be N x 67
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User input (only needed for grayscale image calculations).
% The used values here worked well for the applications in the paper, but
% this may depend on the size of the data clusters.
UpscaleFactor = 10; % Upscale the coordinates to preserve details.
PadFactor = 3; % Number of pixels to pad the grayscale image.
DilationFactor = 1; % Degree to which dilation/erosion should occur.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
if size(Data{1},2) > 2
    Data = cellfun(@(x) x(:,1:2),Data,'UniformOutput',false);
end

% Determine the order of the descriptors to be calculated.
% These are divided into several groups, depending on their method of
% calculation.
GeometricDescriptors = ["NumberLocs","MajorAxis","MinorAxis","Area","FilledArea","ConvexArea","GyrationRadius","EulerNumber","EigenvaluesRatio","Eigenentropy","Density","AreaRatio","AspectRatio","FormRatio","Rectangularity","Circularity","ConvexCircularity","Convexity","Solidity","EquivalentDiameter","FiberLength","FiberWidth","Curl","MajorAxisNormArea","MinorAxisNormArea","EllipseMajorAxis","EllipseMinorAxis","EllipseAspectRatio","EllipseRatio","EllipseEccentricity"];
BoundaryDescriptors = ["FullPerimeter","OuterPerimeter","ConvexPerimeter","ElasticEnergy","BendingEnergy","MeanCurvature","BendingEnergyNormArea"]; 
SkeletonDescriptors = ["SkeletonCloudWidth","SkeletonTotalLength","SkeletonIntersections","SkeletonMeanLength","SkeletonMeanOrientation","SkeletonMeanTortuosity","SkeletonTotalLengthNormArea","SkeletonIntersectionsNormArea"];
TextureDescriptors = ["Contrast","Correlation","Energy","Entropy","Kurtosis","MeanIntensity","MeanNonZeroIntensity","RMSRoughness","Skewness"];
HuMomentDescriptors = "HuMoment"+arrayfun(@string, 1:7);
FractalDescriptors = ["LocalMinkowskiBouligandDim","MinkowskiSausage","HausdorffDim","BoundaryHausdorffDim","SkeletonMinkowskiSausage","SkeletonHausdorffDim"];

% Combine all descriptors.
DescriptorList = [GeometricDescriptors BoundaryDescriptors SkeletonDescriptors TextureDescriptors HuMomentDescriptors FractalDescriptors]; % Put all shape descriptor names together.

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
    % cloud descriptors
    % ---------------------------------------------------------------------
    
    % Perform a Principal Component Analysis to rotate the data (maximizing
    % the length).
    [~,RotatedCluster,Eigenvalues] = pca(Cluster); % Uses the built-in Matlab pca routine.
    x_length = max(RotatedCluster(:,1))-min(RotatedCluster(:,1)); % Determine the length of the cloud in x.
	y_length = max(RotatedCluster(:,2))-min(RotatedCluster(:,2)); % Determine the length of the cloud in y.

    % Calculate the alphaShape.
    % Extract the full boundary and the boundary of the convex hull of the 
    % point cloud.
    % Extract the polyshape of the point cloud and its border coordinates.
    ClusterAlphaShape = alphaShape(RotatedCluster); % Create the alphaShape (uses Delaunay triangulation).
    AlphaOneRegion = criticalAlpha(ClusterAlphaShape,'one-region'); % Calculate the alpha corresponding to a 'one region' alphaShape, to make sure that all (x,y) coordinates are connected.
    ClusterAlphaShape.Alpha = AlphaOneRegion; % Apply it to the alphaShape.

    I = boundary(ClusterAlphaShape.Points,1); % Calculate the boundary points of the alphaShape and retrieve their indices.
    BorderPoints = ClusterAlphaShape.Points(I,:); % Retrieve the actual boundary points.
    PolygonCluster = polyshape(BorderPoints(:,1),BorderPoints(:,2),'Simplify',false); % Create a polygon of the boundary points (which will be completely filled).

    I = boundary(ClusterAlphaShape.Points,0); % Calculate the convex hull of the alphaShape and retrieve the indices.
    ConvHullPoints = ClusterAlphaShape.Points(I,:); % Retrieve the actual boundary points.
    PolygonConvHullCluster = polyshape(ConvHullPoints(:,1),ConvHullPoints(:,2)); % Create a polygon of the convex hull points (which will be completely filled).

    [BorderCoords, ShapePolygon] = FindBorders(ClusterAlphaShape); % Find inner and outer borders of the point cloud (using alphaShape to save time).

    % Calculate the boundary descriptors.
    ReductionIdx = cell2mat(cellfun(@(x) x(1) >= 5, cellfun(@size, BorderCoords,'UniformOutput',false),'UniformOutput',false)); % Selection based on size.
    BorderCoordsReduced = BorderCoords(ReductionIdx); % Do the actual selection.
    [ElasticEnergy, BendingEnergy] = cellfun(@(x) Energy(x), BorderCoordsReduced,'UniformOutput',false); % Calculate the elastic energy and bending energy for each border.
    ElasticEnergy = sum(vertcat(ElasticEnergy{:})); % The elastic energy of the cluster will be the sum of each border.
    BendingEnergy = sum(vertcat(BendingEnergy{:})); % The bending energy of the cluster will be the sum of each border.

    BorderCurvature = Curvature(BorderPoints); % Only calculate the curvature for the outer boundary of the point cloud. 

    % Store the shape descriptors: localization descriptors.
    DescriptorsTable(i,:).("NumberLocs") = size(RotatedCluster,1);
    DescriptorsTable(i,:).("Area") = area(ClusterAlphaShape);
    DescriptorsTable(i,:).("MajorAxis") = max([x_length,y_length]);
    DescriptorsTable(i,:).("MinorAxis") = min([x_length,y_length]);
    DescriptorsTable(i,:).("FilledArea") = area(PolygonCluster);
    DescriptorsTable(i,:).("ConvexArea") = area(PolygonConvHullCluster);
    DescriptorsTable(i,:).("GyrationRadius") = sqrt(mean((RotatedCluster(:,1)-mean(RotatedCluster(:,1))).^2 + (RotatedCluster(:,2)-mean(RotatedCluster(:,2))).^2));
    DescriptorsTable(i,:).("EulerNumber") = ShapePolygon.NumRegions - ShapePolygon.NumHoles;
    DescriptorsTable(i,:).("EigenvaluesRatio") = Eigenvalues(1) / Eigenvalues(2);
    DescriptorsTable(i,:).("Eigenentropy") = -Eigenvalues(1)/sum(Eigenvalues)*log(Eigenvalues(1)/sum(Eigenvalues)) - Eigenvalues(2)/sum(Eigenvalues)*log(Eigenvalues(2)/sum(Eigenvalues));

    % Store the shape descriptors: boundary descriptors.
    DescriptorsTable(i,:).("FullPerimeter") = perimeter(ClusterAlphaShape);
    DescriptorsTable(i,:).("OuterPerimeter") = perimeter(PolygonCluster);
    DescriptorsTable(i,:).("ConvexPerimeter") = perimeter(PolygonConvHullCluster);
    DescriptorsTable(i,:).("ElasticEnergy") = ElasticEnergy;
    DescriptorsTable(i,:).("BendingEnergy") = BendingEnergy;
    DescriptorsTable(i,:).("MeanCurvature") = mean(BorderCurvature,'omitnan');

    % Store the dependent shape descriptors (coordinate-based).
    DescriptorsTable(i,:).("Density") = DescriptorsTable(i,:).("NumberLocs") / DescriptorsTable(i,:).("Area");
    DescriptorsTable(i,:).("AreaRatio") = DescriptorsTable(i,:).("Area") / DescriptorsTable(i,:).("FilledArea");
    DescriptorsTable(i,:).("AspectRatio") = DescriptorsTable(i,:).("MajorAxis") / DescriptorsTable(i,:).("MinorAxis");
    DescriptorsTable(i,:).("FormRatio") = DescriptorsTable(i,:).("Area") / (DescriptorsTable(i,:).("MajorAxis")^2);
    DescriptorsTable(i,:).("Rectangularity") = DescriptorsTable(i,:).("Area") / (DescriptorsTable(i,:).("MajorAxis")*DescriptorsTable(i,:).("MinorAxis"));
    DescriptorsTable(i,:).("Circularity") = 4 * pi * DescriptorsTable(i,:).("FilledArea") / (DescriptorsTable(i,:).("OuterPerimeter")^2);
    DescriptorsTable(i,:).("ConvexCircularity") = 4 * pi*DescriptorsTable(i,:).("ConvexArea") / (DescriptorsTable(i,:).("ConvexPerimeter")^2);
    DescriptorsTable(i,:).("Convexity") = DescriptorsTable(i,:).("OuterPerimeter") / DescriptorsTable(i,:).("ConvexPerimeter");
    DescriptorsTable(i,:).("Solidity") = DescriptorsTable(i,:).("Area") / DescriptorsTable(i,:).("ConvexArea");
    DescriptorsTable(i,:).("EquivalentDiameter") = sqrt(4 * DescriptorsTable(i,:).("FilledArea") / pi);
    DescriptorsTable(i,:).("FiberLength") = (DescriptorsTable(i,:).("FullPerimeter") + realsqrt(abs((DescriptorsTable(i,:).("FullPerimeter"))^2 - 16 * DescriptorsTable(i,:).("Area")))) / 4;
    DescriptorsTable(i,:).("FiberWidth") = DescriptorsTable(i,:).("Area") / DescriptorsTable(i,:).("FiberLength");
    DescriptorsTable(i,:).("Curl") = DescriptorsTable(i,:).("MajorAxis") / DescriptorsTable(i,:).("FiberLength");
    DescriptorsTable(i,:).("MajorAxisNormArea") = DescriptorsTable(i,:).("MajorAxis") / DescriptorsTable(i,:).("Area");
    DescriptorsTable(i,:).("MinorAxisNormArea") = DescriptorsTable(i,:).("MinorAxis") / DescriptorsTable(i,:).("Area");
    DescriptorsTable(i,:).("BendingEnergyNormArea") = DescriptorsTable(i,:).("BendingEnergy") / DescriptorsTable(i,:).("Area");
    
    % ---------------------------------------------------------------------
    % Calculate the descriptors based on a coordinate mask.
    % Ellipse descriptors - Skeleton descriptors - Texture descriptors - 
    % Fractal descriptors - Humoment descriptors - Dependent descriptors
    % ---------------------------------------------------------------------

    % Caclulate the image transforms needed for the next part.
    ImageTransforms = ImageTransformCloud(ShapePolygon); % Calculate the different images needed to calculate their proprties.
    [FractProperties, SkelProperties] = ImageProperties(ImageTransforms); % Calculate the different properties from the previously obtained images.

    % Grayscale image transforms and their properties.
    [~,GrayscaleProperties] = MakeGrayscaleImage(RotatedCluster, UpscaleFactor, PadFactor, DilationFactor);

    % Calculate central moments: scale invariants and rotation invariants.
    % Scale invariants are used in the Ellipse calculations.
    % Rotation invariants are the HuMoment descriptors.
    Moments = MomentCalculation(ImageTransforms.Logical); % Calculates all moments in a separate function. Technically this can also be done using the rotated cluster coordinates (it should theoretically not change)! It was chosen to implement it this way to avoid numerical issues

    % Store the shape descriptors: Ellipse descriptors.
    DescriptorsTable(i,:).("EllipseMajorAxis") = 2*sqrt(2)*sqrt(Moments.n20 + Moments.n02 + sqrt(Moments.M2));
    DescriptorsTable(i,:).("EllipseMinorAxis") = 2*sqrt(2)*sqrt(Moments.n20 + Moments.n02 - sqrt(Moments.M2));
    DescriptorsTable(i,:).("EllipseAspectRatio") = DescriptorsTable(i,:).("EllipseMajorAxis") / DescriptorsTable(i,:).("EllipseMinorAxis");
    DescriptorsTable(i,:).("EllipseRatio") = DescriptorsTable(i,:).("FilledArea") / (pi * DescriptorsTable(i,:).("EllipseMajorAxis") * DescriptorsTable(i,:).("EllipseMinorAxis") / 4);
    DescriptorsTable(i,:).("EllipseEccentricity") = 2 * sqrt((DescriptorsTable(i,:).("EllipseMajorAxis") / 2)^2 - (DescriptorsTable(i,:).("EllipseMinorAxis") / 2)^2) / DescriptorsTable(i,:).("EllipseMajorAxis");

    % Store the shape descriptors: skeleton descriptors.
    DescriptorsTable(i,:).("SkeletonCloudWidth") = SkelProperties.Cloudwidth;
    DescriptorsTable(i,:).("SkeletonTotalLength") = SkelProperties.TotalLength;
    DescriptorsTable(i,:).("SkeletonIntersections") = SkelProperties.Intersections;
    DescriptorsTable(i,:).("SkeletonMeanLength") = SkelProperties.MeanLength;
    DescriptorsTable(i,:).("SkeletonMeanOrientation") = SkelProperties.MeanOrientation;
    DescriptorsTable(i,:).("SkeletonMeanTortuosity") = SkelProperties.MeanTortuosity;

    % Store the shape descriptors: texture descriptors.
    DescriptorsTable(i,:).("Contrast") = GrayscaleProperties.Contrast;
    DescriptorsTable(i,:).("Correlation") = GrayscaleProperties.Correlation;
    DescriptorsTable(i,:).("Energy") = GrayscaleProperties.Energy;
    DescriptorsTable(i,:).("Entropy") = GrayscaleProperties.Entropy;
    DescriptorsTable(i,:).("Kurtosis") = GrayscaleProperties.Kurtosis;
    DescriptorsTable(i,:).("MeanIntensity") = GrayscaleProperties.MeanIntensity;
    DescriptorsTable(i,:).("MeanNonZeroIntensity") = GrayscaleProperties.MeanNonZeroIntensity;
    DescriptorsTable(i,:).("RMSRoughness") = GrayscaleProperties.RMSRoughness;
    DescriptorsTable(i,:).("Skewness") = GrayscaleProperties.Skewness;

    % Store the shape descriptors: fractal descriptors.
    DescriptorsTable(i,:).("LocalMinkowskiBouligandDim") = FractProperties.MBDim;
    DescriptorsTable(i,:).("MinkowskiSausage") = FractProperties.MSausageDim;
    DescriptorsTable(i,:).("HausdorffDim") = FractProperties.HDim;
    DescriptorsTable(i,:).("BoundaryHausdorffDim") = FractProperties.HDim_Boundary;
    DescriptorsTable(i,:).("SkeletonMinkowskiSausage") = FractProperties.MSausageDim_Skeleton;
    DescriptorsTable(i,:).("SkeletonHausdorffDim") = FractProperties.HDim_Skeleton;

    % Store the shape descriptors: 7 HuMoments.
    for j = 1:7
        DescriptorsTable(i,:).("HuMoment"+string(j)) = Moments.("M"+string(j));
    end

    % Store the dependent shape descriptors.
    DescriptorsTable(i,:).("SkeletonTotalLengthNormArea") = DescriptorsTable(i,:).("SkeletonTotalLength") ./  DescriptorsTable(i,:).("Area");
    DescriptorsTable(i,:).("SkeletonIntersectionsNormArea") = DescriptorsTable(i,:).("SkeletonIntersections") ./  DescriptorsTable(i,:).("Area");

end

% Update the waitbar and then close it.
waitbar(1,wb,'Calculating Descriptors. Finished...');
delete(wb);

% Turn warnings back on
warning('on','all');

end
