function [OrderedBranches, NumberBranchPoints] = ExtractBranches3D(SkeletonImage)
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

% Extract the size of the Skeleton image in x, y, and z.
[m, n, o] = size(SkeletonImage);

% Extract the branchpoints and endpoints for further calculation.
BranchPointImage = bwmorph3(SkeletonImage, 'branchpoints');
EndPointImage = bwmorph3(SkeletonImage, 'endpoints');

% Get the coordinates of the branchpoints and endpoints.
BP_Image = find(BranchPointImage(:));
[bp_x,bp_y,bp_z] = ind2sub([m, n, o],BP_Image);
NumberBranchPoints = numel(BP_Image);

EP_Image = find(EndPointImage(:));
[ep_x,ep_y,ep_z] = ind2sub([m, n, o],EP_Image);

% Construct the Skeleton image without the branch points and get the
% coordinates for these branchpoint-less branches. This allows to extract
% each branch individually using regionprops3.
Branches = SkeletonImage & ~BranchPointImage;
if sum(Branches(:)) == 0
    Branches = SkeletonImage;
end
Stats = regionprops3(Branches,"VoxelList");

% Correct the branch lengths by checking if the list of coordinates
% contains an endpoint. If it does, then the branch was connected to 1
% branchpoint. If it does not, then the branch was connected to 2
% branchpoints. Adjust the branches accordingly.
% This ignores any branches that are only 2 voxels long and those 2 voxels 
% were both branchpoints.
Voxels = Stats.VoxelList;
if ~iscell(Voxels)
    Voxels = {Voxels};
end
ContainsEndPoint = cellfun(@(x) sum(ismember([ep_y ep_x ep_z],x,'rows')),Voxels);
BranchWithEndPointAndBranchPoint = Voxels(ContainsEndPoint==1);
BranchWithTwoBranchPoints = Voxels(ContainsEndPoint==0);
BranchWithTwoEndPoints = Voxels(ContainsEndPoint>1);

Branches = cell(size(Voxels,1),1);
% Situation 1: the branch is one that contains an endpoint and the other
% 'endpoint' is a branchpoint in the Skeleton Image.
if ~isempty(BranchWithEndPointAndBranchPoint)
    
    for i = 1:size(BranchWithEndPointAndBranchPoint,1)

        % Indicate the status of a coordinate (1: endpoint; 2: branchpoint;
        % 0: none of the above)
        BranchWithEndPointAndBranchPoint{i}(:,4) = 0;
        Idx_ep = ismember(BranchWithEndPointAndBranchPoint{i}(:,1:3),[ep_y ep_x ep_z],'rows');
        BranchWithEndPointAndBranchPoint{i}(Idx_ep,4) = 1;

        % Make an image of the branch coordinates (in which the second
        % 'endpoint' was removed). 
        emptyImage = zeros(m,n,o);
        for j = 1:size(BranchWithEndPointAndBranchPoint{i},1)
            emptyImage(BranchWithEndPointAndBranchPoint{i}(j,2),BranchWithEndPointAndBranchPoint{i}(j,1),BranchWithEndPointAndBranchPoint{i}(j,3)) = 1;
        end

        % Figure out which coordinates are endpoints in this new image and
        % remove the original endpoints as we are only interested in the
        % coordinates of the branchpoint that is the other endpoint for 
        % this branch.
        % Two situations: one where the branch is long, and the other where
        % the branch is only 2 voxels long (an endpoint and a branchpoint).
        if size(BranchWithEndPointAndBranchPoint{i},1) > 1
            EndPointsBranch = bwmorph3(emptyImage, 'endpoints');
            for j = 1:size(ep_x,1)
                EndPointsBranch(ep_x(j),ep_y(j),ep_z(j)) = 0;
            end
        else
            EndPointsBranch = emptyImage;
        end

        % Extract the coordinates of the remaining endpoint. The real
        % 'endpoint' (i.e., a branchpoint in the original Skeleton Image)
        % will then be a neighbouring pixel. Convolve with a 3 x 3 x 3 cube
        % and find out which branchpoint of the original Skeleton Image is
        % a part of this new image.
        [EP_Branch_x,EP_Branch_y,EP_Branch_z] = ind2sub([m, n, o],find(EndPointsBranch(:)));
        EndPointsBranch = convn(EndPointsBranch,ones(3,3,3),'same'); % Convolve the remaining endpoint to find the real endpoint.
        [epb_x,epb_y,epb_z] = ind2sub([m, n, o],find(EndPointsBranch(:)));
        BP_Idx = find(ismember([bp_y bp_x bp_z],[epb_y epb_x epb_z],'rows'));
    
        % Situation 1a: Only 1 branchpoint of the original Skeleton Image
        % is found.
        if numel(BP_Idx) == 1
            Branches{i} = vertcat(BranchWithEndPointAndBranchPoint{i},[bp_y(BP_Idx) bp_x(BP_Idx) bp_z(BP_Idx) 2]);
        end
    
        % Situation 1b: More than 1 branchpoint of the original Skeleton
        % Image is found. Pick the closest one. If the multiple 
        % branchpoints are equally close, pick the first one as that 
        % decision does not really matter.
        if numel(BP_Idx) > 1
            DistanceToEndpoint = pdist2([EP_Branch_y EP_Branch_x EP_Branch_z],[bp_y(BP_Idx) bp_x(BP_Idx) bp_z(BP_Idx)]);
            [~,Idx] = min(DistanceToEndpoint);
            BP_Idx = BP_Idx(Idx);
            Branches{i} = vertcat(BranchWithEndPointAndBranchPoint{i},[bp_y(BP_Idx) bp_x(BP_Idx) bp_z(BP_Idx) 2]);
        end

    end

end

% Situation 2: The branch does not contain endpoints. Thus, it contains
% two branchpoints of the original Skeleton Image. Use the same strategy as
% in Situation 1, but now it has to be done for each 'endpoint'.
if ~isempty(BranchWithTwoBranchPoints)

    for i = 1:size(BranchWithTwoBranchPoints,1)

        % Indicate the status of a coordinate (1: endpoint; 2: branchpoint;
        % 0: none of the above)
        BranchWithTwoBranchPoints{i}(:,4) = 0;

        % Make an image of the branch coordinates (in which the second
        % 'endpoint' was removed).
        emptyImage = zeros(m,n,o);
        for j = 1:size(BranchWithTwoBranchPoints{i},1)
            emptyImage(BranchWithTwoBranchPoints{i}(j,2),BranchWithTwoBranchPoints{i}(j,1),BranchWithTwoBranchPoints{i}(j,3)) = 1;
        end

        % Figure out which coordinates are endpoints in this new image.
        % Two situations: one where the branch is long, and the other where
        % the branch is short.
        if size(BranchWithTwoBranchPoints{i},1) > 1
            EndPointsBranch = bwmorph3(emptyImage, 'endpoints');
        else
            EndPointsBranch = emptyImage;
        end

        % Extract the coordinates of the endpoint(s). The real
        % 'endpoint' (i.e., a branchpoint in the original Skeleton Image)
        % will then be a neighbouring pixel. Convolve with a 3 x 3 x 3 cube
        % and find out which branchpoint of the original Skeleton Image is
        % a part of this new image. Do this for each 'endpoint'
        % individually.
        [EP_Branch_x,EP_Branch_y,EP_Branch_z] = ind2sub([m, n, o],find(EndPointsBranch(:)));
        for j = 1:size(EP_Branch_x,1)
            EndPointsBranch = zeros(m,n,o);
            EndPointsBranch(EP_Branch_x(j),EP_Branch_y(j),EP_Branch_z(j)) = 1;
            EndPointsBranch = convn(EndPointsBranch,ones(3,3,3),'same'); % Close the image to avoid having disconnected pixels.
            [epb_x,epb_y,epb_z] = ind2sub([m, n, o],find(EndPointsBranch(:)));

            BP_Idx = cell(size(epb_x,1),1);
            for k = 1:size(epb_x,1)
                BP_Idx{k} = find(ismember([bp_y bp_x bp_z],[epb_y(k) epb_x(k) epb_z(k)],'rows'));
            end
            BP_Idx = vertcat(BP_Idx{:});

            % Situation 2a: Only 1 branchpoint is found.
            if numel(BP_Idx) == 1
                BranchWithTwoBranchPoints{i} = vertcat(BranchWithTwoBranchPoints{i},[bp_y(BP_Idx) bp_x(BP_Idx) bp_z(BP_Idx) 2]);
            end

            % Situation 2b: More than 1 branchpoint is found. Pick the closest
            % one. If the multiple branchpoints are equally close, pick the
            % first one as that decision does not matter then.
            % Situation 2c: More than 1 branchpoint is found, and only 1
            % other point is part of the branch. That means that it is a
            % middle point and it should retain the 2 branch points.
            if numel(BP_Idx) > 1
                if size(BranchWithTwoBranchPoints{i},1) ~= 1
                    DistanceToEndpoint = pdist2([EP_Branch_y EP_Branch_x EP_Branch_z],[bp_y(BP_Idx) bp_x(BP_Idx) bp_z(BP_Idx)]);
                    [~,Idx] = min(DistanceToEndpoint);
                    BP_Idx = BP_Idx(Idx);
                    BranchWithTwoBranchPoints{i} = vertcat(BranchWithTwoBranchPoints{i},[bp_y(BP_Idx) bp_x(BP_Idx) bp_z(BP_Idx) 2]);
                else
                    DistanceToEndpoint = pdist2([EP_Branch_y EP_Branch_x EP_Branch_z],[bp_y(BP_Idx) bp_x(BP_Idx) bp_z(BP_Idx)]);
                    [~,Idx] = mink(DistanceToEndpoint,2);
                    BP_Idx = BP_Idx(Idx);
                    BranchWithTwoBranchPoints{i} = vertcat(BranchWithTwoBranchPoints{i},[bp_y(BP_Idx) bp_x(BP_Idx) bp_z(BP_Idx) 2*ones(numel(BP_Idx),1)]);
                end
            end

        end

        % Add the full branch to the variable.
        Branches{size(BranchWithEndPointAndBranchPoint,1)+i} = BranchWithTwoBranchPoints{i};

    end

end

% Situation 3: The branch contains 2 endpoints of the original Skeleton
% Image.
if ~isempty(BranchWithTwoEndPoints)

    for i = 1:size(BranchWithTwoEndPoints,1)

        % Indicate the status of a coordinate (1: endpoint; 2: branchpoint;
        % 0: none of the above)
        BranchWithTwoEndPoints{i}(:,4) = 0;
        Idx_ep = ismember(BranchWithTwoEndPoints{i}(:,1:3),[ep_y ep_x ep_z],'rows');
        BranchWithTwoEndPoints{i}(Idx_ep,4) = 1;

        % Add the full branch to the variable.
        Branches{size(BranchWithEndPointAndBranchPoint,1)+size(BranchWithTwoBranchPoints,1)+i} = BranchWithTwoEndPoints{i};

    end
    
end

% Branches that are smaller than 2 voxels are always going to be endpoints.
BranchLengths = find(cellfun('size',Branches,1)<2);
if ~isempty(BranchLengths)
    for i = 1:numel(BranchLengths)
        Branches{BranchLengths}(4) = 1;
    end
end

% Order the branches from endpoint to endpoint. This will be done by going
% from one of the endpoints (arbitrary which ones) and calculating the
% distance to each point to find the closest one. And then remove that
% previous point from the pool and do the same thing.
if ~isempty(Branches)
    OrderedBranches = cell(numel(Branches),1);
    for i = 1:numel(Branches)
    
        if size(Branches{i},1) > 1
            % Extract the current branch so it is easier to work with, and find a
            % starting point. Remove that point from the search pool as well.
            CurrentBranch = Branches{i};
            EndPoints = find(CurrentBranch(:,4));
            if ~isempty(EndPoints)
                FirstEndpoint = EndPoints(1);
            else
                FirstEndpoint = 1; % This is when there are no endpoints (a 'closed' loop type of structure).
            end            
            Start = CurrentBranch(FirstEndpoint,:);
            CurrentBranch(FirstEndpoint,:) = [];
        
            % Make the ordered branch and add the starting point to the variable.
            OrderedBranches{i} = NaN(size(Branches{i},1),4);
            OrderedBranches{i}(1,:) = Start;
            for j = 1:size(CurrentBranch,1)-1
        
                % Calculate the distance to the other points and find the one that
                % is closest. This can be done because the branches have been
                % extracted already and so no confusion is possible.
                DistanceToRemainingPoints = pdist2(Start(1:3),CurrentBranch(:,1:3));
                [~,Idx] = min(DistanceToRemainingPoints);
        
                % Add that new point to the ordered branch variable and take it as
                % a new starting point while also removing it from the search pool.
                OrderedBranches{i}(1+j,:) = CurrentBranch(Idx,:);
                Start = CurrentBranch(Idx,:);
                CurrentBranch(Idx,:) = [];
            end
            % Finish the branch with the final point that is left.
            OrderedBranches{i}(end,:) = CurrentBranch;
        else
            OrderedBranches{i} = Branches{i};
        end
    end
else
    OrderedBranches = [];
end

end