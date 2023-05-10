function PCAPlots(Data,Class,D2orD3,varargin)
% -------------------------------------------------------------------------
% Function that plots the PCA plots of a data, and colours the different
% data points according to the different classes.
% Examples on how to use it:
%   PCAPlots(X,Y,'2D');
%   PCAPlots(X,[],'2D',AxisView=5);
%   PCAPlots(X,Y,'3D',SaveAs='path\to\save\image.png',FullScreen=0, ...
%       Title='ExampleTitle',Legend={'Class 1','Class 2'}, ...
%       AxisView=[-10 10 -10 10 -10 10],Ax1LabelPos=[10,2,3,5], ...
%       Ax2LabelPos=[-20,5,8,10]);
% Please note that the syntax on how to specify option input has changed
% since Matlab 2021a. Example of before Matlab 2021a:
%   PCAPlots(X,Y,'2D','FullScreen',1,'Title','ExampleTitle');
% -------------------------------------------------------------------------
% Input:
%   Data:           A matrix containing the data that should be plotted in 
%                   the PCA space
%                       Its size should be m x n (m: number of samples; n: 
%                       number of features)
%   Class:          The classes/groups of the different samples. This can 
%                   be an empty matrix as well (no groups in the data)
%                       Its size should be m x 1 (m: number of samples).
%   D2orD3:         2-dimensional or 3-dimensional PCA plot of the data.
%                       It is specified as: '2D' or '3D'
%
% Optional input:
%   AutoScale:      Select whether or not the data has to be autoscaled
%                   before doing the PCA.
%   SaveAs:         The save name. The figure will be saved as .png and 
%                   .fig. If none is specified, the figure will not save.
%   FullScreen:     Fullscreen (1) figure or not (0). Specified as a scalar
%   Title:          The title of the figure. Specified as a char
%   Legend:         The legend of the different data groups. Specified as a
%                   cell of strings/chars
%   MarkerSize:     The size of each point in the plot (default: 25)
%   AxisView:       How the plot should be viewed. This can be specified as
%                   a scalar (e.g., 5 - represents percent) or a vector 
%                   containing the axis range. This does NOT concern the
%                   azimuth and elevation angles.
%                       If the scalar is used, it will calculate outliers
%                       w.r.t. the median and then largen the FOV with the
%                       percent specified.
%   Ax1LabelPos:    A vector containing the rotation and position (x,y,z) 
%                   of the first axis label (x-axis) (only in 3D plots)
%                       Specified as: [Rotation,x_pos,y_pos,z_pos]
%   Ax2LabelPos:    A vector containing the rotation and position (x,y,z) 
%                   of the second axis label (y-axis) (only in 3D plots)
%                       Specified as: [Rotation,x_pos,y_pos,z_pos]
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

% Set default values if not specified in the input.
DefaultAutoscale = 0;
DefaultSave = [];
DefaultFullScreen = 0;
DefaultTitle = [];
DefaultLegend = [];
DefaultMarkerSize = 50;
DefaultAxisView = [];
DefaultAx1Label = [];
DefaultAx2Label = [];

% Parse the input
p = inputParser;
ValidPlotType = @(x) ischar(x) && (strcmp(x,'2D') || strcmp(x,'3D'));
validScalar1 = @(x) isnumeric(x) && isscalar(x) && ~(x < 0) && ~(x > 1);
validScalar2 = @(x) isnumeric(x) && isscalar(x) && x > 0;
validAxis = @(x) (isnumeric(x) && isscalar(x) && x > 0 && x <= 50) || (isnumeric(x) && isvector(x) && (numel(x) == 4 || numel(x) == 6));
addRequired(p,'Data',@ismatrix);
addRequired(p,'Class',@ismatrix);
addRequired(p,'D2orD3',ValidPlotType);
addOptional(p,'Autoscale',DefaultAutoscale,validScalar1);
addOptional(p,'SaveAs',DefaultSave,@ischar);
addOptional(p,'FullScreen',DefaultFullScreen,validScalar1);
addOptional(p,'Title',DefaultTitle,@ischar);
addOptional(p,'Legend',DefaultLegend,@iscell);
addOptional(p,'MarkerSize',DefaultMarkerSize,validScalar2);
addOptional(p,'AxisView',DefaultAxisView,validAxis);
addOptional(p,'Ax1LabelPos',DefaultAx1Label);
addOptional(p,'Ax2LabelPos',DefaultAx2Label);
parse(p,Data,Class,D2orD3,varargin{:});

Data = p.Results.Data;
Class = p.Results.Class;
D2orD3 = p.Results.D2orD3;
Autoscale = p.Results.Autoscale;
SaveAs = p.Results.SaveAs;
FullScreen = p.Results.FullScreen;
Title = p.Results.Title;
Legend = p.Results.Legend;
MarkerSize = p.Results.MarkerSize;
AxisView = p.Results.AxisView;
Ax1LabelPos = p.Results.Ax1LabelPos;
Ax2LabelPos = p.Results.Ax2LabelPos;

clear p DefaultSave DefaultFullScreen DefaultTitle DefaultLegend DefaultMarkerSize DefaultAxis DefaultAx1Label DefaultAx2Label varargin

% Set default colour/shape scheme
Colors = [1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1;0.8 0.6 0.4; 0 0.4 0;0 0 0;0.5 0.5 0.5;0.4 0.6 0.8];
Shape = {'d','s','^','v','p','o','h','+','*','x'};

% Set some of the parameters
if isempty(Class)
    Class = ones(size(Data,1),1);
end

if size(Class,1) == 1
    Class = Class';
end

if isempty(Title)
    DisplayTitle = 0;
else
    DisplayTitle = 1;
end

if isempty(Legend)
    DisplayLegend = 0;
else
    DisplayLegend = 1;
end

if isempty(Ax1LabelPos)
    ChangeViewAx1 = 0;
else
    ChangeViewAx1 = 1;
end

if isempty(Ax2LabelPos)
    ChangeViewAx2 = 0;
else
    ChangeViewAx2 = 1;
end

if isempty(SaveAs)
    SavePlot = 0;
else
    SavePlot = 1;
end

if Autoscale
    Data = AutoScale(Data);
end

% Calculate the PCA
[U,S] = svds(Data,size(Data,2));
Scores = U(:,1:3)*S(1:3,1:3);
VarExplained = diag(S).^2 / sum(diag(S).^2) * 100;

% Set the axis range for the plot
if ~isempty(AxisView) && isscalar(AxisView)
    Scores_OutliersRemoved = rmoutliers(Scores,'median');
    minScores = min(Scores_OutliersRemoved);
    maxScores = max(Scores_OutliersRemoved);

    AxisRange = zeros(2,3);
    for i = 1:3
        if minScores(i) > 0
            AxisRange(1,i) = minScores(i)*(1-AxisView/100);
        else
            AxisRange(1,i) = minScores(i)*(1+AxisView/100);
        end
        if maxScores(i) > 0
            AxisRange(2,i) = maxScores(i)*(1+AxisView/100);
        else
            AxisRange(2,i) = maxScores(i)*(1-AxisView/100);
        end
    end
    
    Ind = [];
    for i = 1:size(Scores)
        if Scores(i,1) < AxisRange(1,1) || Scores(i,1) > AxisRange(2,1) || Scores(i,2) < AxisRange(1,2) || Scores(i,2) > AxisRange(2,2) || Scores(i,3) < AxisRange(1,3) || Scores(i,3) > AxisRange(2,3)
            Ind = [Ind;i];
        end
    end
    Scores(Ind,:) = [];
    Class(Ind) = [];

    AxisRange = AxisRange(:);
elseif ~isempty(AxisView) && numel(AxisView) > 1
    AxisRange = AxisView;
else
    AxisRange = [];
end

% Open the figure
if FullScreen == 1
    figure('units','normalized','outerposition',[0 0 1 1]);
    set(gcf,'color','white','InvertHardCopy','off')
else
    figure;
    set(gcf,'color','white','InvertHardCopy','off')
end

% Plot the data
[~,UniqueClasses] = groupcounts(Class);
for i = 1:numel(UniqueClasses)
    if strcmp(D2orD3,'2D')
        scatter(Scores(Class==UniqueClasses(i),1),Scores(Class==UniqueClasses(i),2),MarkerSize,Colors(i,:),'filled',Shape{i});
        hold on;
    elseif strcmp(D2orD3,'3D')
        scatter3(Scores(Class==UniqueClasses(i),1),Scores(Class==UniqueClasses(i),2),Scores(Class==UniqueClasses(i),3),MarkerSize,Colors(i,:),'filled',Shape{i});
        hold on;
    end
end
set(gca,'FontWeight','bold');

% Adjust the axis and draw the zero lines
if strcmp(D2orD3,'2D')
    if ~isempty(AxisRange)
        axis(AxisRange(1:4));
    end
    Axis = axis;
    line([Axis(1) Axis(2)],[0 0],[0 0],'Color','k','LineStyle','--','LineWidth',2);
    line([0 0],[Axis(3) Axis(4)],[0 0],'Color','k','LineStyle','--','LineWidth',2);
elseif strcmp(D2orD3,'3D')
    if ~isempty(AxisRange)
        axis(AxisRange);
    end
    Axis = axis;
    line([Axis(1) Axis(2)],[0 0],[0 0],'Color','k','LineStyle','--','LineWidth',2);
    line([0 0],[Axis(3) Axis(4)],[0 0],'Color','k','LineStyle','--','LineWidth',2);
    line([0 0],[0 0],[Axis(5) Axis(6)],'Color','k','LineStyle','--','LineWidth',2);
end

% Show the title if needed
if DisplayTitle == 1
    title(Title,'FontSize',12,'FontWeight','bold');
end

% Show the legend if needed
if DisplayLegend == 1
    legend(Legend,'FontSize',12,'FontWeight','bold','location','best');
end

% Set the axis labels
if strcmp(D2orD3,'2D')
    xlabel(['PC 1 (Variance explained: ' num2str(round(VarExplained(1),2)) '%)']);
    ylabel(['PC 2 (Variance explained: ' num2str(round(VarExplained(2),2)) '%)']);
elseif strcmp(D2orD3,'3D')
    if ChangeViewAx1 == 1
        xlabel(['PC 1 (Variance explained: ' num2str(round(VarExplained(1),2)) '%)'],'Rotation',Ax1LabelPos(1),'Position',Ax1LabelPos(2:4));
    else
        xlabel(['PC 1 (Variance explained: ' num2str(round(VarExplained(1),2)) '%)']);
    end
    if ChangeViewAx2 == 1
        ylabel(['PC 2 (Variance explained: ' num2str(round(VarExplained(2),2)) '%)'],'Rotation',Ax2LabelPos(1),'Position',Ax2LabelPos(2:4));
    else
        ylabel(['PC 2 (Variance explained: ' num2str(round(VarExplained(2),2)) '%)']);
    end
    zlabel(['PC 3 (Variance explained: ' num2str(round(VarExplained(3),2)) '%)'])
end

% Save the figure if needed
if SavePlot == 1
    [p,f,~] = fileparts(SaveAs);
    SaveAs = fullfile(p,f);

    print([SaveAs '.png'],'-dpng','-r300');
    savefig([SaveAs '.fig'])
end