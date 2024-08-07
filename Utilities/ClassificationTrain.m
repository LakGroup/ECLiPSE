function ClassificationModel = ClassificationTrain(Data,Class,Method,varargin)
% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------
% Input:
%   Data:   The X-block of the data (the descriptors)
%               Its size will be m x n (m: number of samples; n: variables)
%   Class:  The Y-block of the data (the classes)
%               Its size will be m x 1 (m: number of samples)
%   Method: The variable selection method
%           Choices:    'knn' (k-nearest neighbours classification from the Statistics & Machine Learning Toolbox; The MathWorks Inc.)
%                       'RandomForest' (Random Forest classification from the Statistics & Machine Learning Toolbox; The MathWorks Inc.)
%                       'PLS-DA' (PLS-Discriminant Analysis from the PLS Toolbox; Eigenvector Inc.)
%                       'LReg-DA' (Logistic Regression-Discriminant Analysis from the PLS Toolbox; Eigenvector Inc.)
%                       'AdaBoostM1' (Adaptive Boosting classification for 2-class problems from the Statistics & Machine Learning Toolbox; The MathWorks Inc.)
%                       'AdaBoostM2' (Adaptive Boosting classification for multiclass problems from the Statistics & Machine Learning Toolbox; The MathWorks Inc.)
%                       'RUSBoost' (Random Undersampling Boosting Boosting classification for 2-class problems from the Statistics & Machine Learning Toolbox; The MathWorks Inc.)
%                       'LogitBoost' (Adaptive Logistic Regression classification for 2-class problems from the Statistics & Machine Learning Toolbox; The MathWorks Inc.)
%
% Output:
%   VarSelModel:    A structure containing different variables depending on
%                   the variable selection method used.
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
Defaultknn_Neighbours = 0;
Defaultknn_Distance = 'euclidean';
DefaultRF_Trees = 100;
DefaultOptions = [];
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
NumLearningCycles = p.Results.NumLearningCycles;
RatioToSmallest = p.Results.RatioToSmallest;

% Turn warnings off, to not flood the command window
warning('off','all');

% Check if the input Data is in the correct format.
if ~ismatrix(Data) || size(Data,2) == 1
    error("The input Data should be specified as a matrix (m x n: m = number of samples; n = number of variables).")
end

% Check if the input Class is in the correct format.
if ~ismatrix(Class) || size(Class,2) ~= 1
    error("The input Class should be specified as a vector (m x 1: m = number of samples).")
end

% Check if the specified method is part of the list
MethodCheck = ismember(Method, {'knn','RandomForest','PLS-DA','LReg-DA','AdaBoostM1','AdaBoostM2','RUSBoost','LogitBoost'});
if MethodCheck == 0
    error("The method can only be one of following (see 'help Classification' for more information): 'knn','RandomForest','PLS-DA','LReg-DA','AdaBoostM1','AdaBoostM2,'RUSBoost','LogitBoost'")
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
            Options.display = 'off';
            Options.waitbar = 'off';
            Options.cvi = {'vet' [] 1};
        end

        % Make the classification model using the training data
        ClassificationModel = lregda(Data,Class,Options);

    % AdaBoostM1
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

    % AdaBoostM2
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

    % RUSBoost
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

    % LogitBoost
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

    % No classification. This should not occur as there is an error
    % check at the start of this function
    otherwise

        % Make an empty classification model
        ClassificationModel = [];

end

% Turn warnings back on
warning('on','all');

end