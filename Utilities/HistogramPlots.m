function HistogramPlots(Data,Class,varargin)
% -------------------------------------------------------------------------
% Function that plots the histograms of the individual descriptors,
% according to a class.
% Examples on how to use it:
%   HistogramPlots(Data,Class);
%   HistogramPlots(Data,Class,SaveAs='path\to\save\',Autoscale=1, ...
%       VariableNames=VarNames,Bins=100,Range=0.05,PlotConfig=[4 5], ...
%       FullScreen=1,Legend={'Class 1','Class 2'});
% Please note that the syntax on how to specify option input has changed
% since Matlab 2021a. Example of before Matlab 2021a:
%   HistogramPlots(Data,Class,'FullScreen',1,'PlotConfig',[4 5]);
% -------------------------------------------------------------------------
% Input:
%   Data:           A matrix containing the data for each of the variables.
%                   Its size will be m x n (m: number of samples; n: number
%                   of variables).
%   Class:          A vector containing the class relationship of each
%                   sample. Its size is m x 1 (m: number of samples).
%
% Optional Input:
%   SaveAs:         The path to save the images in. The figure will be 
%                   saved as .png. If none is specified, the figure will 
%                   not save.
%   AutoScale:      Select whether or not the data has to be autoscaled
%                   before doing the PCA.
%   VariableNames:  A cell containing the names of each of the variables.
%                   Its size is m x 1 (m: number of samples). If it is not
%                   provided, no labels will be used for it.
%   Bins:           A scalar with how many bins a histogram should have.
%                   Default: 100.
%   Range:          A scalar that determines the range of the plotted
%                   histogram values. This can be provided if the data 
%                   contains some outliers. Default: 0.05 (meaning that
%                   everything between 0.05 - 0.95 of the max value will be
%                   plotted).
%   PlotConfig:     The number of subplots inside the figure. This should
%                   be specified as a [m n] vector (m: number of rows, 
%                   n: number of columns). If none is speficied, every
%                   cluster will be saved separately.
%   FullScreen:     Fullscreen (1) figure or not (0). Specified as a scalar
%   Legend:         The legend of the different data groups. Specified as a
%                   cell of strings/chars
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
DefaultSave = [];
DefaultAutoscale = 0;
DefaultVariableNames = [];
DefaultBins = 100;
DefaultRange = 0.05;
DefaultPlotconfig = [1 1];
DefaultFullScreen = 0;
DefaultLegend = [];

% Parse the input
p = inputParser;
validScalar1 = @(x) isnumeric(x) && isscalar(x) && (x == 0 || x == 1);
validScalar2 = @(x) isnumeric(x) && isscalar(x) && x > 0;
validScalar3 = @(x) isnumeric(x) && isscalar(x) && ~(x < 0) && ~(x > 1);
validPlotConfig = @(x) isvector(x) && numel(x)==2 && all(x > 0);
addRequired(p,'Data',@ismatrix);
addRequired(p,'Class',@isvector);
addOptional(p,'SaveAs',DefaultSave,@ischar);
addOptional(p,'Autoscale',DefaultAutoscale,validScalar1);
addOptional(p,'VariableNames',DefaultVariableNames,@iscell);
addOptional(p,'Bins',DefaultBins,validScalar2);
addOptional(p,'Range',DefaultRange,validScalar3);
addOptional(p,'PlotConfig',DefaultPlotconfig,validPlotConfig);
addOptional(p,'FullScreen',DefaultFullScreen,validScalar1);
addOptional(p,'Legend',DefaultLegend,@iscell);
parse(p,Data,Class,varargin{:});

Data = p.Results.Data;
Class = p.Results.Class;
SaveAs = p.Results.SaveAs;
Autoscale = p.Results.Autoscale;
VariableNames = p.Results.VariableNames;
Bins = p.Results.Bins;
Range = p.Results.Range;
PlotConfig = p.Results.PlotConfig;
FullScreen = p.Results.FullScreen;
Legend = p.Results.Legend;

clear p varargin validScalar2 validScalar2 validScalar3 validPlotConfig DefaultSave DefaultVariableNames DefaultBins DefaultRange DefaultPlotconfig DefaultFullScreen DefaultLegend

% Set default colour/shape scheme
Colors = [1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1;0.8 0.6 0.4; 0 0.4 0;0 0 0;0.5 0.5 0.5;0.4 0.6 0.8];

% Do some check
if ~isempty(VariableNames) && numel(VariableNames) ~= size(Data,2)
    error('The provided size of the VariableNames cell is not the same as the size of the data');
end

if ~isempty(Legend) && numel(Legend) ~= numel(unique(Class))
    error('The provided size of the Legend cell is not the same as the number of classes');
end

% Set some default parameters
if isempty(SaveAs)
    SavePlot = 0;
else
    SavePlot = 1;
end

% Make sure the bins is a scalar
Bins = round(Bins);

% Autoscale the data if needed
if Autoscale
    Data = AutoScale(Data);
end

% To make sure the legend can be seen
if prod(PlotConfig) > 10
    PlotsPerFigure = prod(PlotConfig)-1;
else
    PlotsPerFigure = prod(PlotConfig);
end

% Do the plotting
for i = 1:ceil(size(Data,2)/PlotsPerFigure)

    % Open the figure
    if FullScreen == 1
        figure('units','normalized','outerposition',[0 0 1 1]);
        set(gcf,'color','white','InvertHardCopy','off')
    else
        figure;
        set(gcf,'color','white','InvertHardCopy','off')
    end

    for j = 1:PlotsPerFigure

        % Determine which ones have to be plotted
        k = (i-1)*PlotsPerFigure+j;

        % Reset the 'break' variable that checks whether or not we are past
        % the number of variables or not.
        Break = 0;

        % Plot the clusters
        if ~(k > size(Data,2))

            % Make the histogram
            SortedData = sort(Data(:,k));
            HistIdx = Data(:,k) <= SortedData(round((1-Range)*size(SortedData,1))) & Data(:,k) >= SortedData(round(Range*size(SortedData,1)));
            DataHist = Data(HistIdx,k);
            ClassHist = Class(HistIdx);
            [~,edges] = histcounts(DataHist,Bins,'normalization','probability');

            if size(ClassHist,2) > 1
                ClassHist = ClassHist';
            end
            [~,Groups] = groupcounts(ClassHist);

            subplot(PlotConfig(1),PlotConfig(2),j);
            for l = 1:numel(unique(ClassHist))
                histogram(DataHist(ClassHist==Groups(l)),edges,'normalization','probability','FaceColor',Colors(l,:));hold on;
            end
            hold off

            % Set the axis
            if ~isempty(VariableNames)
                xlabel(VariableNames{k},'FontSize',10,'FontWeight','bold')
            end
            ylabel('Probability (%)','FontSize',10,'FontWeight','bold')
        else
            Break = 1;
            break
        end

    end

    % Set the legend if needed
    if ~isempty(Legend)
        if Break == 1
            subplot(PlotConfig(1),PlotConfig(2),j);
            for l = 1:numel(unique(ClassHist))
                plot(nan, nan, 'Color', Colors(l,:),'LineWidth',10);hold on
                axis off
            end
        elseif PlotsPerFigure > 10
            subplot(PlotConfig(1),PlotConfig(2),j+1);
            for l = 1:numel(unique(ClassHist))
                plot(nan, nan, 'Color', Colors(l,:),'LineWidth',10);hold on
                axis off
            end
        end
        legend(Legend,'FontSize',16,'FontWeight','bold');
    end

    % Save the figure if needed
    if SavePlot
        print([SaveAs num2str(i) '.png'],'-dpng','-r300');
    end

    % Close the open figure
    close
end