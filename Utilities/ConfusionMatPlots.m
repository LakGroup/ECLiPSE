function [ConfusionMat,SumDiags,BestModelIndex] = ConfusionMatPlots(Predicted,GroundTruth,varargin)
% -------------------------------------------------------------------------
% Function that calculates and optionally plots the confusion matrix of the
% classification results.
% Examples on how to use it:
%   ConfusionMatPlots(Predicted,GroundTruth);
%   [ConfusionMat,SumDiags,BestModelIndex] = ConfusionMatPlots( ...
%       Predicted,GroundTruth,SaveAs='path\to\save\',Plot=1, ...
%       FullScreen=1,Title='TitleOfPlot',Labels={'Label1','Label2'}, ...
%       BestNumModels=10,OnlyBestModels=0);
% Please note that the syntax on how to specify option input has changed
% since Matlab 2021a. Example of before Matlab 2021a:
%   ConfusionMatPlots(Predicted,GroundTruth,'SaveAs','path\to\save\', ...
%       'Plot',1);
% -------------------------------------------------------------------------
% Input:
%   Predicted:      A matrix containing the predicted class associations.
%                   Its size will be m x n (m: number of samples, n: number
%                   of models).
%   GroundTruth:    A matrix containing the ground truth class for each
%                   sample.
%                   Its size will be m x n (m: number of samples, n: number
%                   of models).
%
% Optional input:
%   SaveAs:             The path to save the images in. The figure will be 
%                       saved as .png and as .fig. If none is specified, 
%                       the figure will not save.
%   Plot:               Whether or not the confusion matrices just have to
%                       be calculated (0) or also plotted (1). Default: 1.
%   FullScreen:         Fullscreen (1) figure or not (0). Specified as a 
%                       scalar. Default: 0.
%   Title:              The title of the plot. Provided as a char. Default:
%                       [];
%   Labels:             The labels of the individual classes. Default: [].
%   BestNumModels:      Whether or not the best models (based on the
%                       average prediction accuracy) have to be calculated
%                       and/or plotted. Provided as a scalar. Default: 0 
%                       (no calculation).
%   BestNumModelsIdx:   The indices of the best models. Provided as a
%                       vector. Default: [].
%   OnlyBestModels:     Only plot the best models. Provided as a scalar.
%                       Default: 0.
%
% Output:
%   Confusionmat:   The confusion matrices for each of the models
%                   individually. Its size will be (C x C x n; C: number of
%                   classes, n: number of models).
%   SumDiags:       The average prediction accuracy for each model. Its
%                   size will be (n x 1; n: number of models).
%   BestModelIndex: The indices of the best models. Its size will be (s x
%                   1; s: number of best models provided in
%                   'BestNumModels').
% -------------------------------------------------------------------------
% Code written by:
%   Siewert Hugelier    Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   siewert.hugelier@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   xxx
% -------------------------------------------------------------------------

% Set default values if not specified in the input.
DefaultSave = [];
DefaultPlot = 1;
DefaultFullScreen = 0;
DefaultTitle = [];
DefaultLabels = [];
DefaultBestNumModels = 0;
DefaultBestNumModelsIdx = [];
DefaultOnlyBestModels = 0;

% Parse the input
p = inputParser;
validScalar1 = @(x) isnumeric(x) && isscalar(x) && ~(x < 0) && ~(x > 1);
validScalar2 = @(x) isnumeric(x) && isscalar(x) && x > 0;
addRequired(p,'Predicted',@ismatrix);
addRequired(p,'GroundTruth',@ismatrix);
addOptional(p,'Plot',DefaultPlot,validScalar1);
addOptional(p,'SaveAs',DefaultSave,@ischar);
addOptional(p,'FullScreen',DefaultFullScreen,validScalar1);
addOptional(p,'Title',DefaultTitle,@ischar);
addOptional(p,'Labels',DefaultLabels,@iscell);
addOptional(p,'BestNumModels',DefaultBestNumModels,validScalar2);
addOptional(p,'BestNumModelsIdx',DefaultBestNumModelsIdx,@isvector);
addOptional(p,'OnlyBestModels',DefaultOnlyBestModels,validScalar1);
parse(p,Predicted,GroundTruth,varargin{:});

Predicted = p.Results.Predicted;
GroundTruth = p.Results.GroundTruth;
Plot = p.Results.Plot;
SaveAs = p.Results.SaveAs;
FullScreen = p.Results.FullScreen;
Title = p.Results.Title;
Labels = p.Results.Labels;
BestNumModels = p.Results.BestNumModels;
BestNumModelsIdx = p.Results.BestNumModelsIdx;
OnlyBestModels = p.Results.OnlyBestModels;

clear p DefaultSave DefaultFullScreen DefaultTitle DefaultLabels varargin

% Check if the input Predicted is in the correct format.
if ~ismatrix(Predicted)
    error('The input Predicted should be specified as a matrix (m x n: m = number of samples; n = number of models).')
end

% Check if the input Ground truth is in the correct format.
if ~ismatrix(GroundTruth)
    error('The input Ground truth should be specified as a matrix (m x n: m = number of samples; n = number of models).')
end

% Check if the input Predicted and input Ground truth have the same size
if size(Predicted,1) ~= size(GroundTruth,1) || size(Predicted,2) ~= size(GroundTruth,2)
    error(['The input Predicted and input Ground truth should have the same size. Currently: Predicted (' num2str(size(Predicted,1)) 'x' num2str(size(Predicted,2)) ') - Ground Truth (' num2str(size(GroundTruth,1)) 'x' num2str(size(GroundTruth,2)) ').']);
end

% Check if the the best model index combinations are valid
if ((BestNumModels == 0 && isempty(BestNumModelsIdx)) && OnlyBestModels == 1) || (~isempty(BestNumModelsIdx) && BestNumModels == 1)
    error('You did not specify a valid combination to displaying the best models. Please recheck optional parameters ''BestNumModels'', ''BestNumModelsIdx'', and ''OnlyBestModels''.');
end

% Check if the number of best models asked for is higher than the total
% number of models
if BestNumModels > size(GroundTruth,2)
    error(['The number of best models asked for (' num2str(BestNumModels) ' models) is higher than the the total number of models (' num2str(size(GroundTruth,2)) ' models).']);
end

% Check if the number of best models asked for is higher than the total
% number of models
if numel(BestNumModelsIdx) > size(GroundTruth,2)
    error(['The number of best models asked for (' num2str(numel(BestNumModelsIdx)) ' models) is higher than the the total number of models (' num2str(size(GroundTruth,2)) ' models).']);
end

% Set some of the parameters
if isempty(Title)
    DisplayTitle = 0;
else
    DisplayTitle = 1;
end

if isempty(Labels)
    DisplayLabels = 0;
else
    DisplayLabels = 1;
end

if isempty(SaveAs)
    SavePlot = 0;
else
    SavePlot = 1;
end

% Set default colours
nColours = 256;
Cmap1 = [linspace(1,0,nColours)' linspace(1,0,nColours)' ones(nColours,1)];
Cmap2 = [ones(nColours,1) linspace(1,0,nColours)' linspace(1,0,nColours)'];

% Do some initializations
CountsClasses = groupcounts(GroundTruth);
NumModels = size(Predicted,2);
ConfusionMat = zeros(numel(CountsClasses),numel(CountsClasses),NumModels);
SumDiags = zeros(NumModels,1);

% Do the calculations
for i = 1:NumModels
    ConfusionMat(:,:,i) = confusionmat(GroundTruth(:,i),Predicted(:,i));
    SumDiags(i) = sum(diag(ConfusionMat(:,:,i))) ./ numel(GroundTruth(:,i))*100;
    
    ConfusionMat(:,:,i) = bsxfun(@rdivide, ConfusionMat(:,:,i), sum(ConfusionMat(:,:,i),2))*100;

end

if Plot == 1
    % Average the results
    meanConfusionMat = mean(ConfusionMat,3);
    stdConfusionMat = std(ConfusionMat,[],3);

    % Extract the best x models in case this is asked.
    if BestNumModels ~= 0 && isempty(BestNumModelsIdx)
        [~,BestModelIndex] = sort(SumDiags,'descend');
        BestModelIndex = BestModelIndex(1:BestNumModels);
        meanBestConfusionMat = mean(ConfusionMat(:,:,BestModelIndex),3);
        stdBestConfusionMat = std(ConfusionMat(:,:,BestModelIndex),[],3);
    elseif ~isempty(BestNumModelsIdx)
        BestNumModels = numel(BestNumModelsIdx);
        BestModelIndex = BestNumModelsIdx;
        meanBestConfusionMat = mean(ConfusionMat(:,:,BestNumModelsIdx),3);
        stdBestConfusionMat = std(ConfusionMat(:,:,BestNumModelsIdx),[],3);
    else
        BestModelIndex = [];
    end

    % Open the figure
    if FullScreen == 1
        figure('units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color','white','InvertHardCopy','off')
    else
        figure('units','normalized','outerposition',[0.3 0.25 0.6 0.65]);
        set(gcf,'color','white','InvertHardCopy','off')
    end

    % Do the actual drawing
    if OnlyBestModels ~= 1
        if ~isempty(BestModelIndex)
            subplot(1,2,1)
        end
        for i = 1:numel(CountsClasses)
            for j = 1:numel(CountsClasses)
                % Get the text to display
                if size(Predicted,2) >= 2
                    t = {sprintf("%.02f",meanConfusionMat(i,j)),"\pm",sprintf("%.02f%%",stdConfusionMat(i,j))};
                else
                    t = {sprintf("%.02f%%",meanConfusionMat(i,j))};
                end

                % Get the colour maps for the rectangles
                if i==j
                    c = Cmap1(int64(meanConfusionMat(i,j)/100*(nColours-1))+1,:);
                elseif strcmp(t,"(0.00\newline\pm\newline0.00)%")
                    c = [1 1 1];
                    t = "";
                else
                    c = Cmap2(int64(meanConfusionMat(i,j)/100*(nColours-1))+1,:);
                end

                % Draw the rectangles and add the text to it
                rectangle('Position',[j-1,numel(CountsClasses)-i,1,1],'FaceColor',c);
                if (c(1)*76 + c(2)*150 + c(3)*29) > 186
                    textColour = [0 0 0];
                else
                    textColour = [1 1 1];
                end
                t = text(j-0.5,numel(CountsClasses)-i+0.5,t,'HorizontalAlignment','center','FontSize',20,'FontWeight','bold','Color',textColour,'HorizontalAlignment','center');
                if min(c) < 0.4
                    set(t, 'color', 'w');
                end
            end

            axis off
            axis square

            if DisplayLabels == 1
                LabelDisplayed = split(Labels{i});
                text(-0.1,numel(CountsClasses)-i+0.5,LabelDisplayed,'FontSize',20,'FontWeight','bold','HorizontalAlignment','right');
                text(i-0.5,-0.075,LabelDisplayed,'FontSize',20,'FontWeight','bold','HorizontalAlignment','center');
            end
        end
        if DisplayTitle == 1
            title({Title,['Acc: ' num2str(round(mean(SumDiags),2)) ' ± ' num2str(round(std(SumDiags),2)) ' %']},'FontSize',14,'FontWeight','bold')
        end
    end
    if ~isempty(BestModelIndex)
        if OnlyBestModels ~= 1
            subplot(1,2,2)
        end
        for i = 1:numel(CountsClasses)
            for j = 1:numel(CountsClasses)
                % Get the text to display
                if BestNumModels > 1
                    t = {sprintf("%.02f",meanBestConfusionMat(i,j)),"\pm",sprintf("%.02f%%",stdBestConfusionMat(i,j))};
                else
                    t = {sprintf("%.02f%%",meanBestConfusionMat(i,j))};
                end

                % Get the colour maps for the rectangles
                if i==j
                    c = Cmap1(int64(meanBestConfusionMat(i,j)/100*(nColours-1))+1,:);
                elseif strcmp(t,"(0.00\newline\pm\newline0.00)%")
                    c = [1 1 1];
                    t = "";
                else
                    c = Cmap2(int64(meanBestConfusionMat(i,j)/100*(nColours-1))+1,:);
                end

                % Draw the rectangles and add the text to it
                rectangle('Position',[j-1,numel(CountsClasses)-i,1,1],'FaceColor',c);
                if (c(1)*76 + c(2)*150 + c(3)*29) > 186
                    textColour = [0 0 0];
                else
                    textColour = [1 1 1];
                end
                t = text(j-0.5,numel(CountsClasses)-i+0.5,t,'HorizontalAlignment','center','FontSize',20,'FontWeight','bold','Color',textColour,'HorizontalAlignment','center');
                if min(c) < 0.4
                    set(t, 'color', 'w');
                end
            end

            axis off
            axis square

            if DisplayLabels == 1
                LabelDisplayed = split(Labels{i});
                text(-0.1,numel(CountsClasses)-i+0.5,LabelDisplayed,'FontSize',20,'FontWeight','bold','HorizontalAlignment','right');
                text(i-0.5,-0.075,LabelDisplayed,'FontSize',20,'FontWeight','bold','HorizontalAlignment','center');
            end
        end
        if DisplayTitle == 1
            title({[Title ' - Best ' num2str(BestNumModels) ' models'],['Acc: ' num2str(round(mean(SumDiags(BestModelIndex)),2)) ' ± ' num2str(round(std(SumDiags(BestModelIndex)),2)) ' %']},'FontSize',14,'FontWeight','bold')
        end
    end

    if SavePlot == 1
        [p,f,~] = fileparts(SaveAs);
        SaveAs = fullfile(p,f);

        print([SaveAs '.png'],'-dpng','-r300');
        savefig([SaveAs '.fig'])
    end
else
    if BestNumModels ~= 0 && isempty(BestNumModelsIdx)
        [~,BestModelIndex] = sort(SumDiags,'descend');
        BestModelIndex = BestModelIndex(1:BestNumModels);
    elseif ~isempty(BestNumModelsIdx)
        BestModelIndex = BestNumModelsIdx;
    else
        BestModelIndex = [];
    end
end
end