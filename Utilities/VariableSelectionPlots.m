function VariablesSelected = VariableSelectionPlots(Model,varargin)
% -------------------------------------------------------------------------
% Function that visualizes the results obtained from the automatic variable
% selection routine. The point of this routine is to visualize the selected
% variables obtained from multiple individual models (of the same type).
% There are checks in place to ensure that all models come from the same
% variable selection routine (different methods work differently and they
% cannot be treated the same).
% Examples on how to use it:
%   SelectedVariables = VariableSelectionPlots(Boruta_model)
%   SelectedVariables = VariableSelectionPlots(bipls_model,...
%       Threshold=0.6,VariableNames=VarNames,... 
%       SaveAs='path\to\save\image.png');
% Please note that the syntax on how to specify option input has changed
% since Matlab 2021a. Example of before Matlab 2021a:
%   SelectedVariables = VariableSelectionPlots(bipls_model,'Threshold',...
%       0.6,'VariableNames',VarNames,'SaveAs','path\to\save\image.png');
% -------------------------------------------------------------------------
% Input:
%   Model:          The model of the variable selection routine
%                       This will be a structure containing three fields 
%                       (model for individual run; selected variables for 
%                       indidivual runs; method used)
%
% Optional input:
%   Threshold:      An optional parameter to change the threshold (default:
%                   0.5)
%   VariableNames:  An optional parameter containing the names of the
%                   variables (default: 1:#variables)
%   SaveAs:         The name of how to save the file as.
%
% Output:
%   VariablesSelected:  The selected variables for all runs combined (i.e.,
%                       the ones that were selected x % of the time)
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
DefaultThreshold = 0.5;
DefaultSaveAs = [];
DefaultVariableNames = [];

% Parse the input
p = inputParser;
validScalar1 = @(x) isnumeric(x) && isscalar(x) && ~(x < 0) && ~(x > 1);
addRequired(p,'Model',@iscell);
addOptional(p,'Threshold',DefaultThreshold,validScalar1);
addOptional(p,'SaveAs',DefaultSaveAs,@ischar);
addOptional(p,'VariableNames',DefaultVariableNames,@iscell);

parse(p,Model,Threshold,D2orD3,varargin{:});

Model = p.Results.Model;
Threshold = p.Results.Threshold;
SaveAs = p.Results.SaveAs;
VariableNames = p.Results.VariableNames;

% Set the default option for the variable names if it is not provided.
if isempty(VariableNames)
    VariableNames = num2cell(1:size(Model{1}.VarSel));
end

% Check if the input Data is in the correct format.
if ~iscell(Model) || size(Model,2) ~= 1
    error("The input Model should be specified as a cell (m x 1: m = number of samples).")
end

% Check if all provided individual models are the same ones.
Methods = cellfun(@(x) x.Method,Model,'UniformOutput',false);
UniqueMethods = unique(Methods);
if size(UniqueMethods,1) > 1
    error("You cannot mix and match method types. All provided methods have to be the same.")
end

% Check if the specified method is part of the list
MethodCheck = ismember(UniqueMethods, {'biPLS','rPLS','biPLS_PLSToolbox','rPLS_PLSToolbox','iPLS_PLSToolbox','GA_PLSToolbox','ChiSquare','MRMR','ReliefF','Boruta'});
if MethodCheck == 0
    error("The method can only be one of following (see 'help VariableSelection' for more information): 'biPLS','rPLS','biPLS_PLSToolbox','rPLS_PLSToolbox','iPLS_PLSToolbox','GA_PLSToolbox','ChiSquare','MRMR','ReliefF','Boruta'")
end

% Extract the variables selected from each individual model
VarSel = cell2mat(cellfun(@(x) x.VarSel,Model,'UniformOutput',false));

% Calculate the results
MeanOccurance = mean(VarSel);
VariablesSelected = find(MeanOccurance>=Threshold);
NonSelected = MeanOccurance .* (MeanOccurance<Threshold);
Selected = MeanOccurance .* (MeanOccurance>=Threshold);

% Plot the results
figure('units','normalized','outerposition',[0 0 1 1]);
set(gcf,'color','white','InvertHardCopy','off')
bar(NonSelected,'Facecolor',[0.7 0.7 0.7]);
hold on;
bar(Selected,'g');
set(gca,'FontWeight','Bold','XTick',1:size(VariableNames,2),'XTickLabels',VariableNames,'LineWidth',2);
ylabel('Occurance in optimal model (%)','FontWeight','Bold')
line([0 size(VariableNames,2)+1],[Threshold Threshold],'Color','k','LineWidth',2,'LineStyle','--')
title(UniqueMethods,'FontSize',16,'FontWeight','bold')

% Save the figure if needed
if ~isempty(SaveAs)
    [p,f,~] = fileparts(SaveAs);
    SaveAs = fullfile(p,f);

    print([SaveAs '.png'],'-dpng','-r300');
    savefig([SaveAs '.fig'])
end

end