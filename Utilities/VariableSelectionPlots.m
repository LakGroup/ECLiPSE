function VariablesSelected = VariableSelectionPlots(Model,Threshold,VariableNames,SaveAs)
% -------------------------------------------------------------------------
% Function that visualizes the results obtained from the automatic variable
% selection routine. The point of this routine is to visualize the selected
% variables obtained from multiple individual models (of the same type).
% There are checks in place to ensure that all models come from the same
% variable selection routine (different methods work differently and they
% cannot be treated the same).
% Examples on how to use it:
%   SelectedVariables = VariableSelectionPlots(bipls_model,0.6)
%   SelectedVariables = VariableSelectionPlots(Boruta_model)
% -------------------------------------------------------------------------
% Input:
%   Model:          The model of the variable selection routine
%                       This will be a structure containing three fields 
%                       (model for individual run; selected variables for 
%                       indidivual runs; method used)
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

% Set the default option for the threshold if it is not provided.
if nargin < 2 || isempty(Threshold)
    Threshold = 0.5;
end

% Set the default option for the variable names if it is not provided.
if nargin < 3
    VariableNames = num2cell(1:size(Model{1}.VarSel));
end

if nargin < 4
    SaveAs = [];
end

% Check if the input Data is in the correct format.
if ~iscell(Model) || size(Model,2) ~= 1
    error("The input Model should be specified as a cell (m x 1: m = number of samples).")
end

% Check if the input Class is in the correct format.
if ~isscalar(Threshold) || Threshold < 0 || Threshold > 1
    error("The Threshold should be a scalar between 0-1 (both included).")
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