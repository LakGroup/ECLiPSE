function [Branches, NumberBranchPoints] = ExtractBranches(SkeletonImage)
% -------------------------------------------------------------------------
% Function that extracts the branches of a skeletonized image of an object.
% It starts from pixels that are end points in the image and then follows a
% trajectory to the next end point or branch point. The code is not
% optimized for speed (it can probably be improved).
% There is currently also a small bug in there that appears sometimes (not
% always - <0.001% of the clusters analyzed in the work) when only 2 pixels
% are left in the image. These clusters were just removed from the pool.
% Examples on how to use it:
%   [AllBranches,BP_nmb] = ExtractBranches(Skeleton);
% -------------------------------------------------------------------------
% Input:
%   SkeletonImage:  The skeleton image of a pointcloud (obtained using
%                   bwmorp(LogicalImage,'skel',Inf) or
%                   bwskel(LogicalImage).
%
% Output:
%   Branches:           A cell containing the different branches of the
%                       skeleton.
%   NumberBranchPoints: The number of branchpoints of the skeleton.
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

% Set up the lookup tables. The persistent variables are local to a
% function but global in the sense that they are retained in memory.
% The making of these look up tables is slow, so having them as persistent
% variables is a great way to deal with this.
persistent lutep lutbp

if isempty(lutep)
    lutep = makelut(@Endpoint_fcn, 3); % Obtained this LUT from ISBN: 978-0-9820854-0-0 (page 507 and on). This LUT works better than the built-in 'bwmorph(Im,'endpoints')' function. The function they provide is slow, so the lut itself was just copied here.
end
if isempty(lutbp)
    lutbp = makelut(@(x) (sum(x(:)) == 5 & x(5) == 1 & (all(x([1, 3, 7, 9])) | all(x([2, 4, 6, 8])))) | (sum(x(:)) == 4 & x(5) == 1 & (all(x([1, 3, 7])) | all(x([1, 3, 9])) | all(x([1, 7, 9])) | all(x([3, 7, 9])) | all(x([2, 4, 6])) | all(x([2, 4, 8])) | all(x([2, 6, 8])) | all(x([4, 6, 8])) | all(x([1, 6, 8])) | all(x([3, 4, 8])) | all(x([2, 6, 7])) | all(x([2, 4, 9])) | all(x([1, 3, 8])) | all(x([1, 6, 7])) | all(x([3, 4, 9])) | all(x([2, 7, 9])))), 3); % Obtained this LUT from https://github.com/mathworks/detectBranchpoints (this LUT works better than the built-in 'bwmorph(Im,'branchpoints')' function).
end

% Ensure that the skeleton image is thinned because sometimes branchpoints
% are awkward in the skeletonized image.
SkeletonImage = bwmorph(SkeletonImage,'thin',Inf);

% Initialization of the branches.
Branches = cell(0,1); % Create an empty column cell to add the individual branches to.
NumberBranchPoints = 0;

% Extract the images of the different objects in the skeleton image
% The way the polyshapes of the data clouds in this application are set up
% should only give 1 connected component, but this function will also work
% if there are multiple objects in the same field of view.
Objects = bwconncomp(SkeletonImage); % Extract the information about the different objects in the image.
IndividualImages = cell(size(Objects.PixelIdxList,2),1); % Make an empty cell that has the size of the number of objects (should always just be 1, but just in case).
IndividualImages(:) = {zeros(Objects.ImageSize)}; % Fill each of these cells with zeros to create an empty image with the same size as the original image.
for i = 1:size(Objects.PixelIdxList,2)
    IndividualImages{i}(Objects.PixelIdxList{i}) = 1; % Fill up each cell with an image of the skeleton object.
end

% Loop over the different objects and extract the individual branches.
for i = 1:size(IndividualImages,1)

    % Initialization.
    Im = IndividualImages{i}; % Initialize the first image.
    Counter = 0;

    % Branches will never change in the image, so calculate them just once.
    BranchpointImage = bwlookup(Im,lutbp); % Branchpoint image using the lookup table previously defined.
    [BProw, BPcol] = find(BranchpointImage); % Extract the coordinates of the branchpoints.
    BranchPoints = [BProw, BPcol]; % Make a variable containing all the coordinates (easier to work with).
    BranchPointsIDs = sub2ind(size(Im),BranchPoints(:,1),BranchPoints(:,2)); % Convert the coordinates to image elements (needed for ismember).

    % Count the number of branchpoints.
    NumberBranchPoints = NumberBranchPoints + sum(BranchpointImage(:)); % Add the number of branchpoints to the already existing variable.

    % While there are still pixels active in the skeleton image, extract
    % the branches.
    while sum(Im(:)) ~= 0
    
        % Calculate the endpoints of skeleton image at each new iteration.
        % This is because these can change (e.g., branchpoints can become
        % endpoints). If they change from endpoint to branchpoint (i.e.,
        % this means only one branch is still connected to them), then
        % remove them from the variable that keeps track of the branch 
        % points.
        EndpointImage = bwlookup(Im,lutep); % Endpoints image using the lookup table previously defined.
        [EProw, EPcol] = find(EndpointImage); % Extract the coordinates of the endpoints.
        EndPointsIDs = sub2ind(size(Im),EProw,EPcol); % Convert the coordinates to image elements (needed for ismember).

        if ~isempty(BranchPoints)
            BranchPoints = BranchPoints(~ismember(BranchPointsIDs,EndPointsIDs),:); % Check whether the calculated endpoints have coordinates in common with the branchpoints, and reduce the branchpoint variable if necessary.
            BranchPointsIDs = sub2ind(size(Im),BranchPoints(:,1),BranchPoints(:,2)); % Convert the updated coordinates to image elements (needed for ismember).
        end
        
        % The key coordinates are either endpoints or branchpoints and are
        % coordinates from which a branch can start. First we look at the
        % endpoints because less complicated (only 1 possible direction),
        % and then the branchpoints. This is also the variable used to
        % check whether or not a branch end is reached.
        KeyCoords = vertcat([EProw, EPcol],BranchPoints); % Concatenate the endpoints and branchpoints.

        % If this variable is empty, but the skeleton image still has
        % branches, then it is probably because of circular shapes.
        % Randomly assign the new starting point to the first active pixel
        % it finds.
        if isempty(KeyCoords)
            clear KeyCoords % This is because otherwise this variable will stay empty.
            [KeyCoords(:,1),KeyCoords(:,2)] = find(Im,1); % Find the first active pixel and use that as a starting coordinate.
        end

        % Start doing the actual tracking of the branch. This is done by
        % initializing the branch using a key coordinate, and then adding
        % pixels until another keypoint is reached (or until no more points 
        % are available in its surroundings). This code is complicated 
        % because bwtraceboundary does E-SE-S-SW-etc. in the 8-fold 
        % connectivity, whereas a E-S-W-N-SE-SW-NW-NE direction would be 
        % preferred. This is circumvented in this way.
        Branch = [KeyCoords(1,1),KeyCoords(1,2)]; % Initialize the starting point.

        % Loop until one of the ending criteria is met. This code is slow
        % as point per point is added. This is also because this code is
        % actually to trace boundaries and it will always close the loop
        % back to the starting point.
        while true

            FollowCoords = bwtraceboundary(Im,Branch(end,:),'E',4); % Check E-S-W-N first to connect.

            % Check if pixels are in its E-S-W-N surroundings.
            if size(FollowCoords,1) > 2
                Branch = vertcat(Branch,FollowCoords(2,:)); % Add the next point to the branch if there are any.
            else
                FollowCoords = bwtraceboundary(Im,Branch(end,:),'E',8); % Check for an 8-fold connectivity.

                % Check if pixels are in its E-S-W-N-SE-SW-NW-NE surroundings.
                if size(FollowCoords,1) > 2
                    Branch = vertcat(Branch,FollowCoords(2,:)); % Add the next point to the branch if there are any.
                else
                    break % No more points left (circle or isolated line), so stop the branch.
                end
            end

            Im(Branch(end-1,1),Branch(end-1,2)) = 0; % Remove the starting point from the image.

            % Stop the branch if an endpoint or branchpoint is reached.
            if ismember(Branch(end,:),KeyCoords(2:end,:),'rows')
                break
            end

        end

        % Update the pixels in the image. Do this for two different cases:
        % when branchpoints are still present, and when not.
        if ~isempty(BranchPoints)
            
            % Only remove this pixel in the image if it was not a
            % branchpoint (i.e., an endpoint).
            if ~ismember(Branch(end,:),BranchPoints,'rows')
                Im(Branch(end,1),Branch(end,2)) = 0;
            end

            % Check whether or not the starting point was a branchpoint. If
            % yes, then put it back as itmay be part of other branches as
            % well.
            if ismember(Branch(1,:),BranchPoints,'rows')
                Im(Branch(1,1),Branch(1,2)) = 1;
            end
        else
            Im(Branch(end,1),Branch(end,2)) = 0; % Remove this pixel as it is definitely not a branchpoint if that variable is empty.
        end
        
        % Do some cleanup.
        Im = bwareaopen(Im,2,8); % Remove pixels (in this case it will always be branchpoints) that are not surrounded by any other active pixels.
        ActivePixels = find(Im); % Find the current active pixels in the image.
        if ~isempty(BranchPoints)
            BranchPoints = BranchPoints(ismember(BranchPointsIDs,ActivePixels),:); % Remove the coordinates that were just removed by the bwareaopen action.
            BranchPointsIDs = sub2ind(size(Im),BranchPoints(:,1),BranchPoints(:,2)); % Update the ID matrix.
        end

        % Check for some odd cases where a branch is just 2 pixels long and
        % both are branchpoints.
        % Randomly shuffle the branchpoints around to try and fix it. If
        % this does not work (only problem cases will remain after 
        % shuffling), then just recheck for branchpoints and/or use the
        % non-keypoints method.
        % In case everything works as expected (the 'else' case), then add
        % the branch to the output variable.
        if ~isempty(BranchPoints) && size(Branch,1) == 2 && sum(ismember(Branch,BranchPoints,'rows')) == 2
            if Counter < 10
                BranchPoints = BranchPoints(randperm(size(BranchPoints,1)),:); % Randomly permute the branchpoints.
                BranchPointsIDs = sub2ind(size(Im),BranchPoints(:,1),BranchPoints(:,2)); % Update the ID matrix.
                Counter = Counter + 1; % Add to the counter.
            else
                BranchpointImage = bwlookup(Im,lutbp); % Branchpoint image using the lookup table previously defined.
                [BProw, BPcol] = find(BranchpointImage); % Extract the coordinates of the branchpoints.
                BranchPoints = [BProw, BPcol]; % Make a variable containing all the coordinates (easier to work with).
                BranchPointsIDs = sub2ind(size(Im),BranchPoints(:,1),BranchPoints(:,2)); % Convert the coordinates to image elements (needed for ismember).
            end
        else
            Counter = 0; % Reset the counter (just in case).
            Branches = vertcat(Branches,{Branch}); % Add the currently detected branch to the output variable.
        end

        % Clean up the branch follow variable.
        clear Branch

    end

end

% Remove any duplicates from all the branches. This is normally not needed,
% but just to be safe.
Branches = cellfun(@(x) unique(x,'rows','stable'),Branches,'UniformOutput',false); % Keep the 'stable' mode so the order does not change.

end