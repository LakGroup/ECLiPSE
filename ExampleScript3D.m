%% Section 1: Prepare the script to run
% This section takes 0.04 mins on an i7-12700H 2.30GHz Windows 10 system %

% Clear the workspace and 'start fresh'
clc;close all;clear
rng(159753) % For reproducible results

% Start the timer
t1 = tic;

% Set the correct folder
% A safeguard is built in for people who run the script per section.
fullDir = mfilename('fullpath');
fileName = mfilename();
fullDir = extractBefore(fullDir,fileName);

if contains(fullDir,'LiveEditorEvaluationHelper') || contains(fileName,'LiveEditorEvaluationHelper')
    fullDir = fileparts(matlab.desktop.editor.getActiveFilename);
end
cd(fullDir)

% Add the Utilities folder to the path
addpath(genpath([pwd filesep 'Utilities' filesep]))
addpath(genpath([pwd filesep 'Utilities' filesep 'Utilities3D']))

% Load the data
load([pwd filesep 'Data' filesep 'ExampleData3D.mat']);

% Stop the timer
disp(['Section 1 finished. Time elapsed: ' num2str(toc(t1)/60) ' minutes'])

% Clean up the Workspace
clear fullDir fileName t1

%% Section 2: Calculate the descriptors for each data
% This section takes 30.48 mins on an i7-12700H 2.30GHz Windows 10 system %

% Start the timer
t1 = tic;

% Make the Descriptors folder if it does not exist.
if ~exist([pwd filesep 'Descriptors' filesep], 'dir')
    mkdir([pwd filesep 'Descriptors' filesep])
end

% Calculate the descriptors for each different group of data.
DescriptorsLysosomes = CalcDescriptors3D(Lysosomes);
DescriptorsMitochondria = CalcDescriptors3D(Mitochondria);

VariableNames = DescriptorsLysosomes.Properties.VariableNames;

% Remove the incorrectly calculated entries, and convert the table to a
% matrix.
% Note that there are no incorrectly calculated entries in this example, 
% but this can be used when some clusters are behaving strangely to remove 
% them from the data
[DescriptorsLysosomes,ClassLysosomes] = RemoveNaNInf(DescriptorsLysosomes,Class=ClassLysosomes);
[DescriptorsMitochondria,ClassMitochondria] = RemoveNaNInf(DescriptorsMitochondria,Class=ClassMitochondria);

save([pwd filesep 'Descriptors' filesep 'AllDescriptors3D.mat'],'DescriptorsLysosomes','ClassLysosomes','DescriptorsMitochondria','ClassMitochondria','VariableNames')

% Stop the timer
disp(['Section 2 finished. Time elapsed: ' num2str(toc(t1)/60) ' minutes'])

% Clean up the Workspace
clear t1

%% Section 3: Explore the data (before variable selection)
% This section takes 0.05 mins on an i7-12700H 2.30GHz Windows 10 system %

% Start the timer
t1 = tic;

% Make the PCAPlots folder if it does not exist
if ~exist([pwd filesep 'PCAPlots' filesep], 'dir')
    mkdir([pwd filesep 'PCAPlots' filesep])
end

% Make the full dataset
Data = vertcat(DescriptorsLysosomes,DescriptorsMitochondria);
Class = vertcat(ClassLysosomes,ClassMitochondria);

% 2D PCA plot
Autoscale = 1;
Legend = {'Lysosomes','Mitochondria'};
SaveAs = [pwd filesep 'PCAPlots' filesep 'PCAPlot_3DData_2D_before_VarSel.png'];
PCAPlots(Data,Class,'2D',Autoscale=Autoscale,SaveAs=SaveAs,FullScreen=1,Legend=Legend,AxisView=5);

% 3D PCA plot
SaveAs = [pwd filesep 'PCAPlots' filesep 'PCAPlot_3DData_3D_before_VarSel.png'];
PCAPlots(Data,Class,'3D',Autoscale=Autoscale,SaveAs=SaveAs,FullScreen=1,Legend=Legend,AxisView=5);

% Stop the timer
disp(['Section 3 finished. Time elapsed: ' num2str(toc(t1)/60) ' minutes'])

% Clean up the Workspace
clear t1 Autoscale Legend SaveAs

%% Section 4: Perform the Variable Selection
% This section takes 8.90 mins on an i7-12700H 2.30GHz Windows 10 system %

% Start the timer
t1 = tic;

% Make the VariableSelection folder if it does not exist
if ~exist([pwd filesep 'VariableSelection' filesep], 'dir')
    mkdir([pwd filesep 'VariableSelection' filesep])
end

% -------------------------------------------------------------------------
% Parameters to change if necessary
VarReps = 50; % Will only be used if VarSel is set to 1. Recommend to put it low-ish (25-50).
TrainFraction = 0.5; % The fraction of data points included in the variable selection dataset
VarSelMethod = {'rPLS','ReliefF','MRMR','Boruta'}; % For choices, see 'help VariableSelection' - Only 4 were selected to speed up the example, but all of them could be selected as well.
% -------------------------------------------------------------------------

% Do some initializations
VarSelModel = cell(VarReps,numel(VarSelMethod));
CountClasses = groupcounts(Class);
TrainNumber = floor(min(CountClasses)*TrainFraction);

% Perform the variable selection
wb = waitbar(0,['Performing variable selection. Please be patient... At iteration: 0/' num2str(VarReps) ' '],'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(wb,'canceling',0);
for i = 1:VarReps

    % Decide what to do when cancel is pressed.
    if getappdata(wb,'canceling')
        break
    end

    % Update the waitbar
    waitbar(i/VarReps,wb,['Performing variable selection. Please be patient... At iteration: ' num2str(i) '/' num2str(VarReps)]);

    % Create the training data
    [TrainData,TrainClass] = CreateTrainData(Data,Class,ClassSize=TrainNumber);

    % Perform the variable selection
    for j = 1:numel(VarSelMethod)
        VarSelModel{i,j} = VariableSelection(TrainData,TrainClass,VarSelMethod{j});
    end

    % Clear the work space
    clear TrainData TrainClass

end

% Close the waitbar
delete(wb);

% Save the results
for i = 1:numel(VarSelMethod)
    VarSel = VarSelModel(:,i);
    save([pwd filesep 'VariableSelection' filesep 'VariableSelection_3DData_' VarSelMethod{i} '.mat'],'VarSel');
end
save([pwd filesep 'VariableSelection' filesep 'VariableSelection_3DData_AllModels.mat'],'VarSelModel');

% Stop the timer
disp(['Section 4 finished. Time elapsed: ' num2str(toc(t1)/60) ' minutes'])

% Clean up the Workspace
clear VarReps VarSelPlots CountClasses i j t1 TrainNumber TrainFraction wb VarSel

%% Section 5: Plot the Variable Selection results
% This section takes 0.11 mins on an i7-12700H 2.30GHz Windows 10 system %

% Start the timer
t1 = tic;

% Make the VariableSelection folder if it does not exist
if ~exist([pwd filesep 'VariableSelection' filesep], 'dir')
    mkdir([pwd filesep 'VariableSelection' filesep])
end

% -------------------------------------------------------------------------
% Parameters to change if necessary
Threshold_Individual = [0.5 0.5 0.2 0.5];
Threshold_Global = 2;
% -------------------------------------------------------------------------

% Plot the plots for each of the methods
SelectedVariables = cell(size(VarSelModel,2),1);
for i = 1:size(VarSelModel,2)
    SaveAs = [pwd filesep 'VariableSelection' filesep 'Results_3DData_' VarSelModel{1,i}.Method '.png'];
    SelectedVariables{i} = VariableSelectionPlots(VarSelModel(:,i),Threshold=Threshold_Individual(i),VariableNames=VariableNames,SaveAs=SaveAs);
end

% Extract the variables that were selected each time.
AllVariablesSelected = cell2mat(SelectedVariables');
[CountsVariablesSelected,VariablesSelected] = groupcounts(AllVariablesSelected');

% Prepare the results to be plotted
SelectedVars = zeros(size(VariableNames,2),1);
SelectedVars(VariablesSelected) = CountsVariablesSelected;

% Make the colormap so that the results look nice
clrmap = [repmat([0.7 0.7 0.7],[Threshold_Global-1 1]); repmat([0 1 0],[size(VarSelModel,2)-(Threshold_Global-1) 1])];

% do the actual plotting
figure('units','normalized','outerposition',[0 0 1 1]); set(gcf,'color','white','InvertHardCopy','off')
for i = 1:size(VarSelModel,2)
    bar(SelectedVars.*(SelectedVars==i),'FaceColor',clrmap(i,:))
    hold on
end
set(gca,'FontSize',16,'FontWeight','Bold','XTick',[1 5:5:65],'LineWidth',2);
ylabel('Occurance in all methods','FontSize',24,'FontWeight','Bold')
% Draw a line that indicates which threshold was used.
line([0 size(VariableNames,2)+1],[Threshold_Global-0.5 Threshold_Global-0.5],'Color','k','LineWidth',2,'LineStyle','--')
% Save the plot
print([pwd filesep 'VariableSelection' filesep 'Results_3DData_AllMethods.png'],'-dpng','-r300');
savefig([pwd filesep 'VariableSelection' filesep 'Results_3DData_AllMethods.fig'])

% Perform the final variable selection by only using the ones that were
% selected multiple times
VarSel = find(SelectedVars>=Threshold_Global);

% Stop the timer
disp(['Section 5 finished. Time elapsed: ' num2str(toc(t1)/60) ' minutes'])

% Save the results
save([pwd filesep 'VariableSelection' filesep 'VariablesSelected_3DData.mat'],'VarSel')

% Clean up the Workspace
clear t1 Threshold_Individual Threshold_Global SelectedVariables i AllVariablesSelected CountsVariablesSelected VariablesSelected SelectedVars clrmap VarSelModel SaveAs VarSelMethod

%% Section 6: Explore the data (after variable selection)
% This section takes 0.22 mins on an i7-12700H 2.30GHz Windows 10 system %

% Start the timer
t1 = tic;

% Make the PCAPlots folder if it does not exist
if ~exist([pwd filesep 'PCAPlots' filesep], 'dir')
    mkdir([pwd filesep 'PCAPlots' filesep])
end

% 2D PCA plot
Autoscale = 1;
Legend = {'Lysosomes','Mitochondria'};
SaveAs = [pwd filesep 'PCAPlots' filesep 'PCAPlot_3DData_2D_after_VarSel.png'];
PCAPlots(Data(:,VarSel),Class,'2D',Autoscale=Autoscale,SaveAs=SaveAs,FullScreen=1,Legend=Legend,AxisView=5);

% 3D PCA plot
SaveAs = [pwd filesep 'PCAPlots' filesep 'PCAPlot_3DData_3D_after_VarSel.png'];
PCAPlots(Data(:,VarSel),Class,'3D',Autoscale=Autoscale,SaveAs=SaveAs,FullScreen=1,Legend=Legend,AxisView=5);

% Plot the histogram results of the variable selection
if ~exist([pwd filesep 'PCAPlots' filesep 'VarSel_3DData' filesep 'Selected' filesep], 'dir')
    mkdir([pwd filesep 'PCAPlots' filesep 'VarSel_3DData' filesep 'Selected' filesep])
end
SaveAs = [pwd filesep 'PCAPlots' filesep 'VarSel_3DData' filesep 'Selected' filesep];
Legend = {'Lysosomes','Mitochondria'};
HistogramPlots(Data(:,VarSel),Class,SaveAs=SaveAs,Autoscale=1,VariableNames=VariableNames(VarSel),Bins=100,Range=0.05,PlotConfig=[4 5],FullScreen=1,Legend=Legend);

if ~exist([pwd filesep 'PCAPlots' filesep 'VarSel_3DData' filesep 'NotSelected' filesep], 'dir')
    mkdir([pwd filesep 'PCAPlots' filesep 'VarSel_3DData' filesep 'NotSelected' filesep])
end
SaveAs = [pwd filesep 'PCAPlots' filesep 'VarSel_3DData' filesep 'NotSelected' filesep];
NotVarSel = setdiff(1:size(Data,2),VarSel);
HistogramPlots(Data(:,NotVarSel),Class,SaveAs=SaveAs,Autoscale=1,VariableNames=VariableNames(NotVarSel),Bins=100,Range=0.05,PlotConfig=[4 5],FullScreen=1,Legend=Legend);

% Stop the timer
disp(['Section 6 finished. Time elapsed: ' num2str(toc(t1)/60) ' minutes'])

% Clean up the Workspace
clear t1 Autoscale Legend SaveAs

%% Section 7: Automatic classification
% This section takes 30.68 mins on an i7-12700H 2.30GHz Windows 10 system %

% Start the timer
t1 = tic;

% Make Classification folder if it does not exist
if ~exist([pwd filesep 'Classification' filesep], 'dir')
    mkdir([pwd filesep 'Classification' filesep])
end

% -------------------------------------------------------------------------
% Parameters to change if necessary
ClassReps = 200; % Number of reptitions of the model
ClassMethod = {'knn','RandomForest','LogitBoost'}; % For choices, see 'help ClassificationTrain'
TrainFraction = 0.5; % The fraction of data points included in the training dataset for classification
% -------------------------------------------------------------------------

% Set up some initialization
CountsClasses = groupcounts(Class);
TrainClassSize = floor(min(CountsClasses)*TrainFraction);
ValidationClass = zeros((min(CountsClasses)*numel(CountsClasses))-TrainClassSize*numel(CountsClasses),ClassReps);

Prediction_KNN = zeros((min(CountsClasses)*numel(CountsClasses))-TrainClassSize*numel(CountsClasses),ClassReps);
VS_Prediction_KNN = zeros((min(CountsClasses)*numel(CountsClasses))-TrainClassSize*numel(CountsClasses),ClassReps);
Prediction_RF = zeros((min(CountsClasses)*numel(CountsClasses))-TrainClassSize*numel(CountsClasses),ClassReps);
VS_Prediction_RF = zeros((min(CountsClasses)*numel(CountsClasses))-TrainClassSize*numel(CountsClasses),ClassReps);
Prediction_LogitBoost = zeros((min(CountsClasses)*numel(CountsClasses))-TrainClassSize*numel(CountsClasses),ClassReps);
VS_Prediction_LogitBoost = zeros((min(CountsClasses)*numel(CountsClasses))-TrainClassSize*numel(CountsClasses),ClassReps);

Model_KNN = cell(ClassReps,1);
VS_Model_KNN = cell(ClassReps,1);
Model_LogitBoost = cell(ClassReps,1);
VS_Model_LogitBoost = cell(ClassReps,1);
% Note that the models for Random Forest will not be saved because they are
% extremely big and would require you to have quite a bit of RAM memory. 

% Perform the classification for the supervised method
wb = waitbar(0,['Performing classification. Please be patient... At iteration: 0/' num2str(ClassReps) ' '],'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(wb,'canceling',0);
for i = 1:ClassReps

    % Decide what to do when cancel is pressed.
    if getappdata(wb,'canceling')
        break
    end

    % Update the waitbar
    waitbar(i/ClassReps,wb,['Performing classification. Please be patient... At iteration: ' num2str(i) '/' num2str(ClassReps)]);

    % Create the training data
    [TrainData,TrainClass,ValidationData,ValidationClass(:,i)] = CreateTrainData(Data,Class,ClassSize=TrainClassSize,NeedTestData=1);

    % Perform the classification
    if any(ismember(ClassMethod,'knn'))
        Model_KNN{i} = ClassificationTrain(TrainData,TrainClass,'knn',knn_Neighbours=5,knn_Distance='seuclidean');
        Prediction_KNN(:,i) = ClassificationPrediction(ValidationData,Model_KNN{i},'knn');
        VS_Model_KNN{i} = ClassificationTrain(TrainData(:,VarSel),TrainClass,'knn',knn_Neighbours=5,knn_Distance='seuclidean');
        VS_Prediction_KNN(:,i) = ClassificationPrediction(ValidationData(:,VarSel),VS_Model_KNN{i},'knn');
    end
    if any(ismember(ClassMethod,'RandomForest'))
        Model = ClassificationTrain(TrainData,TrainClass,'RandomForest',RF_Trees=250);
        Prediction_RF(:,i) = ClassificationPrediction(ValidationData,Model,'RandomForest');
        Model = ClassificationTrain(TrainData(:,VarSel),TrainClass,'RandomForest',RF_Trees=250);
        VS_Prediction_RF(:,i) = ClassificationPrediction(ValidationData(:,VarSel),Model,'RandomForest');
    end
    if any(ismember(ClassMethod,'LogitBoost'))
        Model = ClassificationTrain(TrainData,TrainClass,'LogitBoost');
        Prediction_LogitBoost(:,i) = ClassificationPrediction(ValidationData,Model,'LogitBoost');
        Model = ClassificationTrain(TrainData(:,VarSel),TrainClass,'LogitBoost');
        VS_Prediction_LogitBoost(:,i) = ClassificationPrediction(ValidationData(:,VarSel),Model,'LogitBoost');
    end
    
    % Clear the work space
    clear TrainData TrainClass ValidationData Model

end

if any(ismember(ClassMethod,'knn'))
    save([pwd filesep 'Classification' filesep 'ClassificationResults_3DData_knn.mat'],'Prediction_KNN','VS_Prediction_KNN','ValidationClass','Model_KNN','VS_Model_KNN')
end
if any(ismember(ClassMethod,'RandomForest'))
    save([pwd filesep 'Classification' filesep 'ClassificationResults_3DData_RandomForest.mat'],'Prediction_RF','VS_Prediction_RF','ValidationClass')
end
if any(ismember(ClassMethod,'LogitBoost'))
    save([pwd filesep 'Classification' filesep 'ClassificationResults_3DData_LogitBoost.mat'],'Prediction_LogitBoost','VS_Prediction_LogitBoost','ValidationClass','Model_LogitBoost','VS_Model_LogitBoost')
end

% Close the waitbar
delete(wb);

% Stop the timer
disp(['Section 7 finished. Time elapsed: ' num2str(toc(t1)/60) ' minutes'])

% Clean up the Workspace
clear ClassReps TrainFraction CountsClasses TrainClassSize i Model

%% Section 8: Plot the classification results
% This section takes 0.14 mins on an i7-12700H 2.30GHz Windows 10 system %

% Start the timer
t1 = tic;

% Make Classification folder if it does not exist
if ~exist([pwd filesep 'Classification' filesep], 'dir')
    mkdir([pwd filesep 'Classification' filesep])
end

% Set common values for the data
Label = {'Lysosomes','Mitochondria'};
NumberOfModels = 10;

% Show the confusion matrices for the supervised methods
if any(ismember(ClassMethod,'knn'))
    SaveAs = [pwd filesep 'Classification' filesep 'Prediction_3DData_KNN.png'];
    Title = 'K-nearest Neighbours';
    [~,Acc_KNN,Idx_KNN] = ConfusionMatPlots(Prediction_KNN,ValidationClass,SaveAs=SaveAs,FullScreen=1,Title=Title,Labels=Label,BestNumModels=NumberOfModels);
    BestModels_KNN = Model_KNN(Idx_KNN);
    Acc_KNN = Acc_KNN(Idx_KNN);
    Best_Acc{1,1} = {[num2str(round(mean(Acc_KNN),2)) ' ± ' num2str(round(std(Acc_KNN),2)) '%']};

    SaveAs = [pwd filesep 'Classification' filesep 'VS_Prediction_3DData_KNN.png'];
    Title = 'K-nearest Neighbours - Variable Selection';
    [~,Acc_VS_KNN,Idx_VS_KNN] = ConfusionMatPlots(VS_Prediction_KNN,ValidationClass,SaveAs=SaveAs,FullScreen=1,Title=Title,Labels=Label,BestNumModels=NumberOfModels);
    BestModels_VS_KNN = VS_Model_KNN(Idx_VS_KNN);
    Acc_VS_KNN = Acc_VS_KNN(Idx_VS_KNN);
    Best_Acc{1,2} = {[num2str(round(mean(Acc_VS_KNN),2)) ' ± ' num2str(round(std(Acc_VS_KNN),2)) '%']};
end
if any(ismember(ClassMethod,'RandomForest'))
    SaveAs = [pwd filesep 'Classification' filesep 'Prediction_3DData_RF.png'];
    Title = 'Random Forest (250 Trees)';
    [~,Acc_RF,Idx_RF] = ConfusionMatPlots(Prediction_RF,ValidationClass,SaveAs=SaveAs,FullScreen=1,Title=Title,Labels=Label,BestNumModels=NumberOfModels);
    Acc_RF = Acc_RF(Idx_RF);
    Best_Acc{2,1} = {[num2str(round(mean(Acc_RF),2)) ' ± ' num2str(round(std(Acc_RF),2)) '%']};

    SaveAs = [pwd filesep 'Classification' filesep 'VS_Prediction_3DData_RF.png'];
    Title = 'Random Forest (250 Trees) - Variable Selection';
    [~,Acc_VS_RF,Idx_VS_RF] = ConfusionMatPlots(VS_Prediction_RF,ValidationClass,SaveAs=SaveAs,FullScreen=1,Title=Title,Labels=Label,BestNumModels=NumberOfModels);
    Acc_VS_RF = Acc_VS_RF(Idx_VS_RF);
    Best_Acc{2,2} = {[num2str(round(mean(Acc_VS_RF),2)) ' ± ' num2str(round(std(Acc_VS_RF),2)) '%']};
end
if any(ismember(ClassMethod,'LogitBoost'))
    SaveAs = [pwd filesep 'Classification' filesep 'Prediction_3DData_LogitBoost.png'];
    Title = 'LogitBoost';
    [~,Acc_LogitBoost,Idx_LogitBoost] = ConfusionMatPlots(Prediction_LogitBoost,ValidationClass,SaveAs=SaveAs,FullScreen=1,Title=Title,Labels=Label,BestNumModels=NumberOfModels);
    Acc_LogitBoost = Acc_LogitBoost(Idx_LogitBoost);
    Best_Acc{3,1} = {[num2str(round(mean(Acc_LogitBoost),2)) ' ± ' num2str(round(std(Acc_LogitBoost),2)) '%']};

    SaveAs = [pwd filesep 'Classification' filesep 'VS_Prediction_3DData_LogitBoost.png'];
    Title = 'LogitBoost - Variable Selection';
    [~,Acc_VS_LogitBoost,Idx_VS_LogitBoost] = ConfusionMatPlots(VS_Prediction_LogitBoost,ValidationClass,SaveAs=SaveAs,FullScreen=1,Title=Title,Labels=Label,BestNumModels=NumberOfModels);
    Acc_VS_LogitBoost = Acc_VS_LogitBoost(Idx_VS_LogitBoost);
    Best_Acc{3,2} = {[num2str(round(mean(Acc_VS_LogitBoost),2)) ' ± ' num2str(round(std(Acc_VS_LogitBoost),2)) '%']};
end

% Change the best results to a table format and set the row/column names
Best_Acc = cell2table(Best_Acc);
Best_Acc.Properties.VariableNames = {'Without variable selection','With variable selection'};
Best_Acc.Properties.RowNames = {'K-nearest Neighbours','Random Forest','LogitBoost'};


% Show a summary of the results
disp('A summary of the classification results is found below (average prediction accuracy):')
disp(Best_Acc)

% Stop the timer
disp(['Section 8 finished. Time elapsed: ' num2str(toc(t1)/60) ' minutes'])