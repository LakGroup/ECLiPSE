function ClassificationModel = ClassificationTrain(Data,Class,Method,varargin)
% -------------------------------------------------------------------------
% Function that trains the classification models based on data training
% that is provided. The Class vector associated to it is different for
% supervised (vector of class associations per sample) and unsupervised (a
% scalar representing the number of classes) methods.
% Examples on how to use it:
%   ClassificationModel = ClassificationTrain(Data,Class,'knn');
%   ClassificationModel = ClassificationTrain(Data,Class,'LReg-DA', ...
%       Options=Options_LRegDA);
%   ClassificationModel = ClassificationTrain(Data,Class,'Kmeans', ...
%       DoPCA=1,nComp=5);
% Please note that the syntax on how to specify option input has changed
% since Matlab 2021a. Example of before Matlab 2021a:
%   ClassificationModel = ClassificationTrain(Data,Class,'Kmeans',...
%       'DoPCA',1,'nComp',5);
% -------------------------------------------------------------------------
% Input:
%   Data:   The X-block of the data (the descriptors)
%               Its size will be m x n (m: number of samples; n: variables)
%   Class:  The Y-block of the data (the classes)
%               Its size should be m x 1 (m: number of samples) for
%               supervised classification methods
%               Its size should be 1 x 1 for unsupervised classification
%               methods (i.e., the number of groups in the data)
%   Method: The variable selection method
%           Choices:    'knn' (k-nearest neighbours classification from the Statistics & Machine Learning Toolbox; The MathWorks Inc.)
%                       'RandomForest' (Random Forest classification from the Statistics & Machine Learning Toolbox; The MathWorks Inc.)
%                       'PLS-DA' (PLS-Discriminant Analysis from the PLS Toolbox; Eigenvector Inc.)
%                       'LReg-DA' (Logistic Regression-Discriminant Analysis from the PLS Toolbox; Eigenvector Inc.)
%                       'AdaBoostM1' (Adaptive Boosting classification for 2-class problems from the Statistics & Machine Learning Toolbox; The MathWorks Inc.)
%                       'AdaBoostM2' (Adaptive Boosting classification for multiclass problems from the Statistics & Machine Learning Toolbox; The MathWorks Inc.)
%                       'RUSBoost' (Random Undersampling Boosting Boosting classification for 2-class problems from the Statistics & Machine Learning Toolbox; The MathWorks Inc.)
%                       'LogitBoost' (Adaptive Logistic Regression classification for 2-class problems from the Statistics & Machine Learning Toolbox; The MathWorks Inc.)
%                       'HCA' (Hierarchical Cluster analysis for unsupervised classification from the PLS Toolbox; Eigenvector Inc.)
%                       'Kmeans' (k-means clustering for unsupervised classification from the Statistics & Machine Learning Toolbox; The MathWorks Inc.)
%
% Optional input:
%   knn_Neighbours: The number of neighbours considered in the knn model.
%                   See 'help fitcknn' for more information.
%   knn_Distance:   The distance metric considered in the knn model. See
%                   'help fitcknn' for more information on the available
%                   options.
%   RF_Trees:       The number of trees used in the Random Forest model.
%                   See 'help TreeBagger' for more information.
%   Options:        Options for 'PLS-DA', 'LReg-DA', and 'HCA' that can be 
%                   provided. See 'help lregda' or 'help plsda' for more 
%                   information.
%   DoPCA:          Do PCA before using KMeans or not (1/0 default: 1)
%   nComp:          The number of components to do the PCA with (only 
%                   Kmeans). Default: 5.
%   Learner:        The learner provided for the 'AdaBoostM1',
%                   'AdaBoostM2', 'RUSBoost', and 'LogitBoost' 
%                   classification methods. See 'help fitcensemble' for
%                   more information.
%   NumLearningCycles:  The number of cycles used in the learning process
%                       for 'AdaBoostM1', 'AdaBoostM2', 'RUSBoost', and 
%                       'LogitBoost' classification methods. See 'help
%                       fitcensemble' for more information. Dfault: 1000.e
%   RatioToSmallest:    Sampling portion for 'RUSBoost'. See 'help
%                       fitcensemble' for more information. Default: 1 for
%                       each class.
%
% Output:
%   ClassificationModel:    A structure containing the classification
%                           model. The contents will depend on the
%                           classification method.
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
Defaultknn_Neighbours = 0;
Defaultknn_Distance = 'seuclidean';
DefaultRF_Trees = 100;
DefaultOptions = [];
DefaultDoPCA = 1;
DefaultnComp = 5;
DefaultLearner = [];
DefaultNumLearningCycles = 0;
DefaultRatioToSmallest = [];

% Parse the input
p = inputParser;
validScalar1 = @(x) isnumeric(x) && isscalar(x) && x > 2;
addRequired(p,'Data',@ismatrix);
addRequired(p,'Class',@ismatrix);
addRequired(p,'Method',@ischar);
addOptional(p,'knn_Neighbours',Defaultknn_Neighbours,validScalar1);
addOptional(p,'knn_Distance',Defaultknn_Distance,@ischar);
addOptional(p,'RF_Trees',DefaultRF_Trees,validScalar1);
addOptional(p,'Options',DefaultOptions);
addOptional(p,'DoPCA',DefaultDoPCA);
addOptional(p,'nComp',DefaultnComp);
addOptional(p,'Learner',DefaultLearner);
addOptional(p,'NumLearningCycles',DefaultNumLearningCycles);
addOptional(p,'RatioToSmallest',DefaultRatioToSmallest);
parse(p,Data,Class,Method,varargin{:});

Data = p.Results.Data;
Class = p.Results.Class;
Method = p.Results.Method;
knn_Neighbours = p.Results.knn_Neighbours;
knn_Distance = p.Results.knn_Distance;
RF_Trees = p.Results.RF_Trees;
Options = p.Results.Options;
Learner = p.Results.Learner;
DoPCA = p.Results.DoPCA;
nComp = p.Results.nComp;
NumLearningCycles = p.Results.NumLearningCycles;
RatioToSmallest = p.Results.RatioToSmallest;

% Update the path to make sure the PLS toolbox is in the appropriate
% position, depending on the method being used
if ismember(Method,{'PLS-DA','LReg-DA','HCA'})
    oldPath = path;
    seperatePaths = split(oldPath,pathsep);
    plsOrNot = contains(seperatePaths,'PLS_Toolbox');

    if isempty(plsOrNot)
        error("You need to install the PLS toolbox (Eigenvector Inc.) to be able to do this. You can download a free trial from www.eigenvector.com");
    end

    newPath = join(join(seperatePaths(plsOrNot),pathsep), join(seperatePaths(~plsOrNot),pathsep));
    addpath(newPath{1});
else
    oldPath = path;
    seperatePaths = split(oldPath,pathsep);
    plsOrNot = contains(seperatePaths,'PLS_Toolbox');
    newPath = join(join(seperatePaths(~plsOrNot),pathsep), join(seperatePaths(plsOrNot),pathsep));
    addpath(newPath{1});
end

% Check if the Class vector is empty or not. Make a vector of 1s if none is
% provided. Downstream error checks may still be valid.
if isempty(Class)
    Class = ones(size(Data,1),1);
end

% Turn warnings off, to not flood the command window
warning('off','all');

% Check if the input Data is in the correct format.
if ~ismatrix(Data) || size(Data,2) == 1
    error("The input Data should be specified as a matrix (m x n: m = number of samples; n = number of variables).")
end

% Check if the input Class is in the correct format.
if ismember(Method,{'knn','RandomForest','PLS-DA','LReg-DA','AdaBoostM1','AdaBoostM2','RUSBoost','LogitBoost'}) && (~ismatrix(Class) || size(Class,2) ~= 1)
    error("The input Class should be specified as a vector for supervised classification (m x 1: m = number of samples).")
elseif ismember(Method,{'HCA','Kmeans'}) && (~ismatrix(Class) || size(Class,2) ~= 1)
    error("The input Class should be specified as a vector for supervised classification (m x 1: m = number of samples).")
end

% Check if the specified method is part of the list
if ~ismember(Method, {'knn','RandomForest','PLS-DA','LReg-DA','AdaBoostM1','AdaBoostM2','RUSBoost','LogitBoost','HCA','Kmeans'})
    error("The method can only be one of following (see 'help Classification' for more information): 'knn', 'RandomForest', 'PLS-DA', 'LReg-DA', 'AdaBoostM1', 'AdaBoostM2, 'RUSBoost', 'LogitBoost', 'HCA', 'Kmeans'")
end

% Check if the method is compatible with the data
if numel(unique(Class)) == 1 && ismember(Method,{'knn','RandomForest','PLS-DA','LReg-DA','AdaBoostM2','RUSBoost','AdaBoostM1','LogitBoost'})
    error("The Class vector only contains a single class. No classification can be performed.")
elseif numel(unique(Class)) == 2 && ~ismember(Method,{'AdaBoostM1','LogitBoost'})
    error("The data is a 2-class classification problem. Only 'AdaBoostM1' and 'LogitBoost' can be used as a method for this data.")
elseif numel(unique(Class)) > 2 && ~ismember(Method,{'knn','RandomForest','PLS-DA','LReg-DA','AdaBoostM2','RUSBoost'})
    error("The data is a multiclass classification problem. Only 'knn', 'RandomForest', 'PLS-DA', 'LReg-DA', 'AdaBoostM2', and 'RUSBoost' can be used as a method for this data.")
elseif numel(Class) == 1 && ~ismember(Method,{'HCA','Kmeans'})
    error("You provided a class structure only compatible with unsupervised classification. Only 'HCA' and 'Kmeans' can be used as a method for this data.")
end

% Train the classification model
switch Method

    % K-nearest neighbours classifier. This is a function that is built-in 
    % in the Statistics and Machine Learning toolbox of Matlab.
    case 'knn'

        % Autoscale the data, as it is not done within the knn routine
        Data = AutoScale(Data);

        % Make the classification model using the training data
        if knn_Neighbours == 0
            ClassificationModel = fitcknn(Data,Class,'Distance',knn_Distance);
        else
            ClassificationModel = fitcknn(Data,Class,'Distance',knn_Distance,'NumNeighbors',knn_Neighbours);
        end

    % Random Forest classifier. This is a function that is built-in 
    % in the Statistics and Machine Learning toolbox of Matlab.
    case 'RandomForest'

        % Autoscale the data, as it is not done within the Random Forest
        % routine
        Data = AutoScale(Data);

        % Make the classification model using the training data
%         ClassificationModel = TreeBagger(RF_Trees,Data,Class,'Method','classification','OOBPredictorImportance','On'); % When you want information on the variable selection
        ClassificationModel = TreeBagger(RF_Trees,Data,Class,'Method','classification');

    % Partial Least Squares Discriminant Analysis from the PLS toolbox. A 
    % 30-day trial demo can be downloaded from www.eigenvector.com
    case 'PLS-DA'

        % Create the options if none are specified
        if isempty(Options)
            Options = plsda('options');
            Options.preprocessing = {'autoscale' 'autoscale'};
            Options.plots = 'none';
            Options.display = 'off';
        end

        % Do an initial model with many components
        model = plsda(Data, Class, size(Data,2)-1, Options);
        model = model.crossvalidate(Data);

        % Do the final model with an optimized number of components.
        [~,ncomp] = min(mean(model.detail.classerrcv));
        ClassificationModel = plsda(Data, Class, ncomp, Options);

    % Logistic Regression Discriminant Analysis from the PLS toolbox. A 
    % 30-day trial demo can be downloaded from www.eigenvector.com
    case 'LReg-DA'

        % Create the options if none are specified
        if isempty(Options)
            Options = lregda('options');
            Options.preprocessing = {preprocess('default','autoscale') preprocess('default','meancenter')};
            Options.plots = 'none';
            Options.display = 'off';
            Options.waitbar = 'off';
            Options.compression = 'pls';
            Options.compressncomp = 5;
            Options.cvi = {'vet' 5 1};
        end

        % Make the classification model using the training data
        ClassificationModel = lregda(Data,Class,Options);

    % AdaBoostM1. This is a function that is built-in 
    % in the Statistics and Machine Learning toolbox of Matlab.
    case 'AdaBoostM1'

        % Autoscale the data, as it is not done within the AdaboostM1
        % routine
        Data = AutoScale(Data);

        % Create the options if none are specified
        if isempty(Learner)
            Learner = templateTree('MaxNumSplits',10);
        end

        % Create the options if none are specified
        if NumLearningCycles == 0
            NumLearningCycles = 1000;
        end

        % Make the classification model using the training data
        ClassificationModel = fitcensemble(Data,Class,'Method','AdaBoostM1','NumLearningCycles',NumLearningCycles,'Learners',Learner,'LearnRate',0.1,'nprint',NumLearningCycles+1);

    % AdaBoostM2. This is a function that is built-in 
    % in the Statistics and Machine Learning toolbox of Matlab.
    case 'AdaBoostM2'

        % Autoscale the data, as it is not done within the AdaBoostM2
        % routine
        Data = AutoScale(Data);

        % Create the options if none are specified
        if isempty(Learner)
            Learner = templateTree('MaxNumSplits',10);
        end

        % Create the options if none are specified
        if NumLearningCycles == 0
            NumLearningCycles = 1000;
        end

        % Make the classification model using the training data
        ClassificationModel = fitcensemble(Data,Class,'Method','AdaBoostM2','NumLearningCycles',NumLearningCycles,'Learners',Learner,'LearnRate',0.1,'nprint',NumLearningCycles+1);

    % RUSBoost. This is a function that is built-in 
    % in the Statistics and Machine Learning toolbox of Matlab.
    case 'RUSBoost'

        % Autoscale the data, as it is not done within the RUSBoost
        % routine
        Data = AutoScale(Data);

        % Create the options if none are specified
        if isempty(Learner)
            Learner = templateTree('MaxNumSplits',10);
        end

        % Create the options if none are specified
        if NumLearningCycles == 0
            NumLearningCycles = 1000;
        end

        % Create the options if none are specified
        if isempty(RatioToSmallest)
            RatioToSmallest = ones(numel(groupcounts(Class)),1);
        end

        % Make the classification model using the training data
        ClassificationModel = fitcensemble(Data,Class,'Method','RUSBoost','NumLearningCycles',NumLearningCycles,'Learners',Learner,'LearnRate',0.1,'RatioToSmallest',RatioToSmallest,'nprint',NumLearningCycles+1);

    % LogitBoost. This is a function that is built-in 
    % in the Statistics and Machine Learning toolbox of Matlab.
    case 'LogitBoost'

        % Autoscale the data, as it is not done within the LogitBoost
        % routine
        Data = AutoScale(Data);

        % Create the options if none are specified
        if isempty(Learner)
            Learner = templateTree('MaxNumSplits',10);
        end

        % Create the options if none are specified
        if NumLearningCycles == 0
            NumLearningCycles = 1000;
        end

        % Make the classification model using the training data
        ClassificationModel = fitcensemble(Data,Class,'Method','LogitBoost','NumLearningCycles',NumLearningCycles,'Learners',Learner,'LearnRate',0.1,'nprint',NumLearningCycles+1);

    % Hierarchical Cluster Analysis from the PLS toolbox (unsupervised). A 
    % 30-day trial demo can be downloaded from www.eigenvector.com
    case 'HCA'

        % Create the options if none are specified
        if isempty(Options)
            Options = cluster('options');
            Options.preprocessing = {'autoscale' 'autoscale'};
            Options.plots = 'none';
            Options.distance = 'euclidean';
            Options.ncomp = size(Data,2)/2;
            Options.pca = 'on';
            Options.mahalanobis = 'off';
        end

        % Make the classification model using the training data
        ClassificationModel = cluster(Data, Options);
        Threshold = mean(ClassificationModel.dist(end-(Class-1):end-(Class-2)));
        ClassificationModel = ClassificationModel.class(find(ClassificationModel.dist<Threshold,1,'last'),:);

        % Order the class groups from 1:n (n being the number of classes).
        % This method is not ideal, and will most likely not correspond to
        % the classes initially provided, but this is hard to do logically.
        [~,Groups] = groupcounts(ClassificationModel');
        for i = 1:Class
            ClassificationModel(ClassificationModel == Groups(i)) = i;
        end

    % K-means classification. This is a function that is built-in 
    % in the Statistics and Machine Learning toolbox of Matlab. 
    case 'Kmeans'
        
        % Autoscale the data, as it is not done within the LogitBoost
        % routine
        Data = AutoScale(Data);

        % Perform PCA if needed
        if DoPCA
            [U,S] = svds(Data,size(Data,2));
            Data = U(:,1:nComp)*S(1:nComp,1:nComp);
        end

        % Create the options
        Dist = 'sqeuclidean';
        Reps = 5;

        % Make the classification model using the training data
        [~,ClassificationModel] = kmeans(Data,Class,'Distance',Dist,'Replicates',Reps);

    % No classification. This should not occur as there is an error
    % check at the start of this function
    otherwise

        % Make an empty classification model
        ClassificationModel = [];

end

% Turn warnings back on
warning('on','all');

end