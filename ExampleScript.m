%% Section 1: Prepare the script to run
% This section takes 0.01 mins on an i7-12700H 2.30GHz Windows 10 system %

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

% Load the data
load([pwd filesep 'Data' filesep 'ExampleData.mat']);

% Stop the timer
disp(['Section 1 finished. Time elapsed: ' num2str(toc(t1)/60) ' minutes'])

% Clean up the Workspace
clear fullDir fileName t1

%% Section 2: Calculate the descriptors for each data
% This section takes 27.5 mins on an i7-12700H 2.30GHz Windows 10 system %

% Start the timer
t1 = tic;

% Make the Descriptors folder if it does not exist.
if ~exist([pwd filesep 'Descriptors' filesep], 'dir')
    mkdir([pwd filesep 'Descriptors' filesep])
end

% Calculate the descriptors for each different group of data.
DescriptorsTau = CalcDescriptors(Tau);
DescriptorsNPCs = CalcDescriptors(NPCs);
DescriptorsMicrotubules = CalcDescriptors(Microtubules);
DescriptorsLysosomes = CalcDescriptors(Lysosomes);
DescriptorsMitochondria = CalcDescriptors(Mitochondria);

VariableNames = DescriptorsTau.Properties.VariableNames;

% Remove the incorrectly calculated entries, and convert the table to a
% matrix.
% Note that there are no incorrectly calculated entries in this example, 
% but this can be used when some clusters are behaving strangely to remove 
% them from the data
[DescriptorsTau,ClassTau] = RemoveNaNInf(DescriptorsTau,Class=ClassTau);
[DescriptorsNPCs,ClassNPCs] = RemoveNaNInf(DescriptorsNPCs,Class=ClassNPCs);
[DescriptorsMicrotubules,ClassMicrotubules] = RemoveNaNInf(DescriptorsMicrotubules,Class=ClassMicrotubules);
[DescriptorsLysosomes,ClassLysosomes] = RemoveNaNInf(DescriptorsLysosomes,Class=ClassLysosomes);
[DescriptorsMitochondria,ClassMitochondria] = RemoveNaNInf(DescriptorsMitochondria,Class=ClassMitochondria);

save([pwd filesep 'Descriptors' filesep 'AllDescriptors.mat'],'DescriptorsTau','ClassTau','DescriptorsNPCs','ClassNPCs','DescriptorsMicrotubules','ClassMicrotubules','DescriptorsLysosomes','ClassLysosomes','DescriptorsMitochondria','ClassMitochondria','VariableNames')

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
Data = vertcat(DescriptorsTau,DescriptorsNPCs,DescriptorsMicrotubules,DescriptorsLysosomes,DescriptorsMitochondria);
Class = horzcat(ClassTau,ClassNPCs,ClassMicrotubules,ClassLysosomes,ClassMitochondria)';

% 2D PCA plot
Autoscale = 1;
Legend = {'Tau','NPCs','Microtubules','Lysosomes','Mitochondria'};
SaveAs = [pwd filesep 'PCAPlots' filesep 'PCAPlot_2D_before_VarSel.png'];
PCAPlots(Data,Class,'2D',Autoscale=Autoscale,SaveAs=SaveAs,FullScreen=1,Legend=Legend,AxisView=5);

% 3D PCA plot
SaveAs = [pwd filesep 'PCAPlots' filesep 'PCAPlot_3D_before_VarSel.png'];
PCAPlots(Data,Class,'3D',Autoscale=Autoscale,SaveAs=SaveAs,FullScreen=1,Legend=Legend,AxisView=5);

% Stop the timer
disp(['Section 3 finished. Time elapsed: ' num2str(toc(t1)/60) ' minutes'])

% Clean up the Workspace
clear t1 Autoscale Legend SaveAs

%% Section 4: Perform the Variable Selection
% This section takes 6.16 mins on an i7-12700H 2.30GHz Windows 10 system %

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
VarSelMethod = {'rPLS','ReliefF','rPLS_PLSToolbox'}; % For choices, see 'help VariableSelection' - Only 2 were selected to speed up the example, but all of them could be selected as well.
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
    save([pwd filesep 'VariableSelection' filesep 'VariableSelection_' VarSelMethod{i} '.mat'],'VarSel');
end
save([pwd filesep 'VariableSelection' filesep 'VariableSelection_AllModels.mat'],'VarSelModel');

% Stop the timer
disp(['Section 4 finished. Time elapsed: ' num2str(toc(t1)/60) ' minutes'])

% Clean up the Workspace
clear VarReps VarSelPlots CountClasses i j t1 TrainNumber TrainFraction wb VarSel

%% Section 5: Plot the Variable Selection results
% This section takes 0.08 mins on an i7-12700H 2.30GHz Windows 10 system %

% Start the timer
t1 = tic;

% Make the VariableSelection folder if it does not exist
if ~exist([pwd filesep 'VariableSelection' filesep], 'dir')
    mkdir([pwd filesep 'VariableSelection' filesep])
end

% -------------------------------------------------------------------------
% Parameters to change if necessary
Threshold_Individual = [0.5 0.5 0.4];
Threshold_Global = 2;
% -------------------------------------------------------------------------

% Plot the plots for each of the methods
SelectedVariables = cell(size(VarSelModel,2),1);
for i = 1:size(VarSelModel,2)
    SaveAs = [pwd filesep 'VariableSelection' filesep 'Results_' VarSelModel{1,i}.Method '.png'];
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
print([pwd filesep 'VariableSelection' filesep 'Results_AllMethods.png'],'-dpng','-r300');
savefig([pwd filesep 'VariableSelection' filesep 'Results_AllMethods.fig'])

% Perform the final variable selection by only using the ones that were
% selected multiple times
VarSel = find(SelectedVars>=Threshold_Global);

% Stop the timer
disp(['Section 5 finished. Time elapsed: ' num2str(toc(t1)/60) ' minutes'])

% Save the results
save([pwd filesep 'VariableSelection' filesep 'VariablesSelected.mat'],'VarSel')

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
Legend = {'Tau','NPCs','Microtubules','Lysosomes','Mitochondria'};
SaveAs = [pwd filesep 'PCAPlots' filesep 'PCAPlot_2D_after_VarSel.png'];
PCAPlots(Data(:,VarSel),Class,'2D',Autoscale=Autoscale,SaveAs=SaveAs,FullScreen=1,Legend=Legend,AxisView=5);

% 3D PCA plot
SaveAs = [pwd filesep 'PCAPlots' filesep 'PCAPlot_3D_after_VarSel.png'];
PCAPlots(Data(:,VarSel),Class,'3D',Autoscale=Autoscale,SaveAs=SaveAs,FullScreen=1,Legend=Legend,AxisView=5);

% Plot the histogram results of the variable selection
if ~exist([pwd filesep 'PCAPlots' filesep 'VarSel' filesep 'Selected' filesep], 'dir')
    mkdir([pwd filesep 'PCAPlots' filesep 'VarSel' filesep 'Selected' filesep])
end
SaveAs = [pwd filesep 'PCAPlots' filesep 'VarSel' filesep 'Selected' filesep];
Legend = {'Tau','NPCs','Microtubules','Lysosomes','Mitochondria'};
HistogramPlots(Data(:,VarSel),Class,SaveAs=SaveAs,Autoscale=1,VariableNames=VariableNames(VarSel),Bins=100,Range=0.05,PlotConfig=[4 5],FullScreen=1,Legend=Legend);

if ~exist([pwd filesep 'PCAPlots' filesep 'VarSel' filesep 'NotSelected' filesep], 'dir')
    mkdir([pwd filesep 'PCAPlots' filesep 'VarSel' filesep 'NotSelected' filesep])
end
SaveAs = [pwd filesep 'PCAPlots' filesep 'VarSel' filesep 'NotSelected' filesep];
NotVarSel = setdiff(1:size(Data,2),VarSel);
HistogramPlots(Data(:,NotVarSel),Class,SaveAs=SaveAs,Autoscale=1,VariableNames=VariableNames(NotVarSel),Bins=100,Range=0.05,PlotConfig=[4 5],FullScreen=1,Legend=Legend);

% Stop the timer
disp(['Section 6 finished. Time elapsed: ' num2str(toc(t1)/60) ' minutes'])

% Clean up the Workspace
clear t1 Autoscale Legend SaveAs

%% Section 7: Automatic classification
% This section takes 12.02 mins on an i7-12700H 2.30GHz Windows 10 system %

% Start the timer
t1 = tic;

% Make Classification folder if it does not exist
if ~exist([pwd filesep 'Classification' filesep], 'dir')
    mkdir([pwd filesep 'Classification' filesep])
end

% -------------------------------------------------------------------------
% Parameters to change if necessary
ClassReps = 200; % Number of reptitions of the model
ClassMethod = {'knn','RandomForest','LReg-DA','Kmeans'}; % For choices, see 'help ClassificationTrain'
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
Prediction_LRegDA = zeros((min(CountsClasses)*numel(CountsClasses))-TrainClassSize*numel(CountsClasses),ClassReps);
VS_Prediction_LRegDA = zeros((min(CountsClasses)*numel(CountsClasses))-TrainClassSize*numel(CountsClasses),ClassReps);

Model_KNN = cell(ClassReps,1);
VS_Model_KNN = cell(ClassReps,1);
Model_LRegDA = cell(ClassReps,1);
VS_Model_LRegDA = cell(ClassReps,1);
Model_Kmeans = cell(1,1);
VS_Model_Kmeans = cell(1,1);
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
    if any(ismember(ClassMethod,'LReg-DA'))
        Model_LRegDA{i} = ClassificationTrain(TrainData,TrainClass,'LReg-DA');
        Prediction_LRegDA(:,i) = ClassificationPrediction(ValidationData,Model_LRegDA{i},'LReg-DA');
        VS_Model_LRegDA{i} = ClassificationTrain(TrainData(:,VarSel),TrainClass,'LReg-DA');
        VS_Prediction_LRegDA(:,i) = ClassificationPrediction(ValidationData(:,VarSel),VS_Model_LRegDA{i},'LReg-DA');
    end
    
    % Clear the work space
    clear TrainData TrainClass ValidationData Model

end

% Perform the unsupervised classification
% Only has to be performed on the data once because the ground truth cannot
% technically be used in this unsupervised process
if any(ismember(ClassMethod,'Kmeans'))
    Model_Kmeans = ClassificationTrain(Data,numel(unique(Class)),'Kmeans');
    Prediction_Kmeans = ClassificationPrediction(Data,Model_Kmeans,'Kmeans');
    VS_Model_Kmeans = ClassificationTrain(Data(:,VarSel),numel(unique(Class)),'Kmeans');
    VS_Prediction_Kmeans = ClassificationPrediction(Data(:,VarSel),VS_Model_Kmeans,'Kmeans');
end

if any(ismember(ClassMethod,'knn'))
    save([pwd filesep 'Classification' filesep 'ClassificationResults_knn.mat'],'Prediction_KNN','VS_Prediction_KNN','ValidationClass','Model_KNN','VS_Model_KNN')
end
if any(ismember(ClassMethod,'RandomForest'))
    save([pwd filesep 'Classification' filesep 'ClassificationResults_RandomForest.mat'],'Prediction_RF','VS_Prediction_RF','ValidationClass')
end
if any(ismember(ClassMethod,'LReg-DA'))
    save([pwd filesep 'Classification' filesep 'ClassificationResults_LRegDA.mat'],'Prediction_LRegDA','VS_Prediction_LRegDA','ValidationClass')%,'Model_LRegDA','VS_Model_LRegDA')
end
if any(ismember(ClassMethod,'Kmeans'))
    save([pwd filesep 'Classification' filesep 'ClassificationResults_Kmeans.mat'],'Prediction_Kmeans','VS_Prediction_Kmeans','ValidationClass','Model_Kmeans','VS_Model_Kmeans')
end

% Close the waitbar
delete(wb);

% Stop the timer
disp(['Section 7 finished. Time elapsed: ' num2str(toc(t1)/60) ' minutes'])

% Clean up the Workspace
clear ClassReps TrainFraction CountsClasses TrainClassSize i Model

%% Section 8: Plot the classification results
% This section takes 6.3 mins on an i7-12700H 2.30GHz Windows 10 system %

% Start the timer
t1 = tic;

% Make Classification folder if it does not exist
if ~exist([pwd filesep 'Classification' filesep], 'dir')
    mkdir([pwd filesep 'Classification' filesep])
end

% Set common values for the data
Label = {'Tau','NPCs','Microtubules','Lysosomes','Mitochondria'};
NumberOfModels = 10;

% Show the confusion matrices for the supervised methods
if any(ismember(ClassMethod,'knn'))
    SaveAs = [pwd filesep 'Classification' filesep 'Prediction_KNN.png'];
    Title = 'K-nearest Neighbours';
    [~,Acc_KNN,Idx_KNN] = ConfusionMatPlots(Prediction_KNN,ValidationClass,SaveAs=SaveAs,FullScreen=1,Title=Title,Labels=Label,BestNumModels=NumberOfModels);
    BestModels_KNN = Model_KNN(Idx_KNN);
    Acc_KNN = Acc_KNN(Idx_KNN);
    Best_Acc{1,1} = {[num2str(round(mean(Acc_KNN),2)) ' ± ' num2str(round(std(Acc_KNN),2)) '%']};

    SaveAs = [pwd filesep 'Classification' filesep 'VS_Prediction_KNN.png'];
    Title = 'K-nearest Neighbours - Variable Selection';
    [~,Acc_VS_KNN,Idx_VS_KNN] = ConfusionMatPlots(VS_Prediction_KNN,ValidationClass,SaveAs=SaveAs,FullScreen=1,Title=Title,Labels=Label,BestNumModels=NumberOfModels);
    BestModels_VS_KNN = VS_Model_KNN(Idx_VS_KNN);
    Acc_VS_KNN = Acc_VS_KNN(Idx_VS_KNN);
    Best_Acc{1,2} = {[num2str(round(mean(Acc_VS_KNN),2)) ' ± ' num2str(round(std(Acc_VS_KNN),2)) '%']};
end
if any(ismember(ClassMethod,'RandomForest'))
    SaveAs = [pwd filesep 'Classification' filesep 'Prediction_RF.png'];
    Title = 'Random Forest (250 Trees)';
    [~,Acc_RF,Idx_RF] = ConfusionMatPlots(Prediction_RF,ValidationClass,SaveAs=SaveAs,FullScreen=1,Title=Title,Labels=Label,BestNumModels=NumberOfModels);
    Acc_RF = Acc_RF(Idx_RF);
    Best_Acc{2,1} = {[num2str(round(mean(Acc_RF),2)) ' ± ' num2str(round(std(Acc_RF),2)) '%']};

    SaveAs = [pwd filesep 'Classification' filesep 'VS_Prediction_RF.png'];
    Title = 'Random Forest (250 Trees) - Variable Selection';
    [~,Acc_VS_RF,Idx_VS_RF] = ConfusionMatPlots(VS_Prediction_RF,ValidationClass,SaveAs=SaveAs,FullScreen=1,Title=Title,Labels=Label,BestNumModels=NumberOfModels);
    Acc_VS_RF = Acc_VS_RF(Idx_VS_RF);
    Best_Acc{2,2} = {[num2str(round(mean(Acc_VS_RF),2)) ' ± ' num2str(round(std(Acc_VS_RF),2)) '%']};
end
if any(ismember(ClassMethod, 'LReg-DA'))
    SaveAs = [pwd filesep 'Classification' filesep 'Prediction_LRegDA.png'];
    Title = 'Logistic Regression';
    [~,Acc_LRegDA,Idx_LRegDA] = ConfusionMatPlots(Prediction_LRegDA,ValidationClass,SaveAs=SaveAs,FullScreen=1,Title=Title,Labels=Label,BestNumModels=NumberOfModels);
    BestModels_LRegDA = Model_KNN(Idx_LRegDA);
    Acc_LRegDA = Acc_LRegDA(Idx_LRegDA);
    Best_Acc{3,1} = {[num2str(round(mean(Acc_LRegDA),2)) ' ± ' num2str(round(std(Acc_LRegDA),2)) '%']};

    SaveAs = [pwd filesep 'Classification' filesep 'VS_Prediction_LRegDA.png'];
    Title = 'Logistic Regression - Variable Selection';
    [~,Acc_VS_LRegDA,Idx_VS_LRegDA] = ConfusionMatPlots(VS_Prediction_LRegDA,ValidationClass,SaveAs=SaveAs,FullScreen=1,Title=Title,Labels=Label,BestNumModels=NumberOfModels);
    BestModels_VS_LRegDA = VS_Model_KNN(Idx_VS_LRegDA);
    Acc_VS_LRegDA = Acc_VS_LRegDA(Idx_VS_LRegDA);
    Best_Acc{3,2} = {[num2str(round(mean(Acc_VS_LRegDA),2)) ' ± ' num2str(round(std(Acc_VS_LRegDA),2)) '%']};
end

% Change the best results to a table format and set the row/column names
Best_Acc = cell2table(Best_Acc);
Best_Acc.Properties.VariableNames = {'Without variable selection','With variable selection'};
Best_Acc.Properties.RowNames = {'K-nearest Neighbours','Random Forest','Logistic Regression'};


% Show a summary of the results
disp('A summary of the classification results is found below (average prediction accuracy):')
disp(Best_Acc)

% Save the images for the unsupervised method
if any(ismember(ClassMethod, 'Kmeans'))

    % Make the K-means folder if it does not exist
    if ~exist([pwd filesep 'Classification' filesep 'Results_Kmeans' filesep], 'dir')
        mkdir([pwd filesep 'Classification' filesep 'Results_Kmeans' filesep])
    end

    % Write some message in the command window
    disp('The results of the K-means unsupervised classification are now saved as images.')

    % Do some initializations
    AllClusters = vertcat(Tau,NPCs,Microtubules,Lysosomes,Mitochondria);
    [~,Groups] = groupcounts(Prediction_Kmeans');
    [~,Groups_VS] = groupcounts(VS_Prediction_Kmeans');
    PlotSizeAll = [16 12 16 10 10];
    PlotSizeAll_VS = [16 12 5 10 10];

    for i = 1:numel(Groups)

        % Make the type folder if it does not exist
        if ~exist([pwd filesep 'Classification' filesep 'Results_Kmeans' filesep 'Type' num2str(i) filesep], 'dir')
            mkdir([pwd filesep 'Classification' filesep 'Results_Kmeans' filesep 'Type' num2str(i) filesep])
        end
        
        % Get the indices of ones that were classified as this type
        CurrentClass = find(Prediction_Kmeans==Groups(i));
        
        % Set the rest of the parameters
        SaveAs = [pwd filesep 'Classification' filesep 'Results_Kmeans' filesep 'Type' num2str(i) filesep];
        PlotSize = PlotSizeAll(i);

        % Do the plotting
        ClusterPlots(AllClusters(CurrentClass),SaveAs=SaveAs,PlotConfig=[5 5],PixelSize=117,PlotSize=PlotSize,FullScreen=1,PlotScalebar=500);
        
    end

    for i = 1:numel(Groups_VS)

        % Make the type folder if it does not exist
        if ~exist([pwd filesep 'Classification' filesep 'Results_Kmeans_VS' filesep 'Type' num2str(i) filesep], 'dir')
            mkdir([pwd filesep 'Classification' filesep 'Results_Kmeans_VS' filesep 'Type' num2str(i) filesep])
        end
        
        % Get the indices of ones that were classified as this type
        CurrentClass = find(VS_Prediction_Kmeans==Groups_VS(i));
        
        % Set the rest of the parameters
        SaveAs = [pwd filesep 'Classification' filesep 'Results_Kmeans_VS' filesep 'Type' num2str(i) filesep];
        PlotSize = PlotSizeAll_VS(i);

        % Do the plotting
        ClusterPlots(AllClusters(CurrentClass),SaveAs=SaveAs,PlotConfig=[5 5],PixelSize=117,PlotSize=PlotSize,FullScreen=1,PlotScalebar=500);
        
    end

end

% Stop the timer
disp(['Section 8 finished. Time elapsed: ' num2str(toc(t1)/60) ' minutes'])