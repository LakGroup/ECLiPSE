function ClusterPlots(Clusters,varargin)
% -------------------------------------------------------------------------
% Function that plots the clusters.
% Examples on how to use it:
%   ClusterPlots(Clusters);
%   ClusterPlots(Clusters,SaveAs='path\to\save\',PlotConfig=[4 5], ...
%       PixelSize=117,PlotSize=10, FullScreen=1, PlotScalebar=500);
% Please note that the syntax on how to specify option input has changed
% since Matlab 2021a. Example of before Matlab 2021a:
%   ClusterPlots(Clusters,'FullScreen',1,'PlotConfig',[4 5]);
% -------------------------------------------------------------------------
% Input:
%   Clusters:       A cell containing the data that should be plotted with
%                   size n x 1
%                       Inside each cell, there should be a matrix with at
%                       least 2 columns (first column: x, second column: y)
%   SaveAs:         The path to save the images in. The figure will be 
%                   saved as .png. If none is specified, the figure will 
%                   not save.
%   PlotConfig:     The number of subplots inside the figure. This should
%                   be specified as a [m n] vector (m: number of rows, 
%                   n: number of columns). If none is speficied, every
%                   cluster will be saved separately.
%   PixelSize:      The size of a single pixel in nm (default: 117)
%   PlotSize:       The size of each plot, specified as a scalar. Default:
%                   10.
%                       10 means that the axis will be: [-10 10 -10 10].
%   FullScreen:     Fullscreen (1) figure or not (0). Specified as a scalar
%   PlotScalebar:   Specify when a scalebar should be plotted. Default: []
%                   The size should be specified in nm. Doublecheck the
%                   PixelSize in this case!
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
DefaultPlotconfig = [1 1];
DefaultPixelSize = 117;
DefaultPlotSize = 10;
DefaultFullScreen = 0;
DefaultPlotScalebar = [];

% Parse the input
p = inputParser;
validScalar1 = @(x) isnumeric(x) && isscalar(x) && x > 0;
validScalar2 = @(x) isnumeric(x) && isscalar(x) && (x == 0 || x == 1);
validPlotConfig = @(x) isvector(x) && numel(x)==2 && all(x > 0);
addRequired(p,'Clusters',@iscell);
addOptional(p,'SaveAs',DefaultSave,@ischar);
addOptional(p,'PlotConfig',DefaultPlotconfig,validPlotConfig);
addOptional(p,'PixelSize',DefaultPixelSize,validScalar1);
addOptional(p,'PlotSize',DefaultPlotSize,validScalar1);
addOptional(p,'FullScreen',DefaultFullScreen,validScalar2);
addOptional(p,'PlotScalebar',DefaultPlotScalebar,validScalar1);
parse(p,Clusters,varargin{:});

Clusters = p.Results.Clusters;
SaveAs = p.Results.SaveAs;
PlotConfig = p.Results.PlotConfig;
PixelSize = p.Results.PixelSize;
PlotSize = p.Results.PlotSize;
FullScreen = p.Results.FullScreen;
PlotScalebar = p.Results.PlotScalebar;

clear p DefaultPlotconfig DefaultSave DefaultPixelSize DefaultPlotSize DefaultPlotScalebar DefaultFullScreen varargin validScalar1 validScalar2 validPlotConfig

% Set some default parameters
if isempty(SaveAs)
    SavePlot = 0;
else
    SavePlot = 1;
end

if isempty(PlotScalebar)
    Scalebar = 0;
else
    Scalebar = PlotScalebar;
end

% Do the plotting
for i = 1:ceil(numel(Clusters)/prod(PlotConfig))

    % Open the figure
    if FullScreen == 1
        figure('units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color','white','InvertHardCopy','off')
    else
        figure;
        set(gcf,'color','white','InvertHardCopy','off')
    end

    for j = 1:prod(PlotConfig)

        % Determine which ones have to be plotted
        k = (i-1)*prod(PlotConfig)+j;

        % Plot the clusters
        if ~(k > numel(Clusters))
            subplot(PlotConfig(1),PlotConfig(2),j);
            plot(Clusters{k}(:,1)-mean(Clusters{k}(:,1)),Clusters{k}(:,2)-mean(Clusters{k}(:,2)),'.k')
    
            % Set the correct axis
            axis([-PlotSize PlotSize -PlotSize PlotSize]);
            axis square
            axis off
    
            % See if a scalebar has to be added or not
            if Scalebar > 0
                line([PlotSize-1-Scalebar/PixelSize PlotSize-1],[-PlotSize+1 -PlotSize+1],'Color','k','LineWidth',5)
                text((PlotSize-1-Scalebar/PixelSize + PlotSize-1)/2,-PlotSize+2,[num2str(Scalebar) ' nm'],'Color','k','FontSize',10,'FontWeight','bold','HorizontalAlignment','center');
            end
        end

    end

    % Save the figure if needed
    if SavePlot
        print([SaveAs num2str(i) '.png'],'-dpng','-r300');
    end

    % Close the open figure
    close
end