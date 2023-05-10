function Prediction = ClassificationPrediction(Data,Model,Method,varargin)
% -------------------------------------------------------------------------
% Function that predicts the class association of each sample based on an
% already existing model. The data should be provided in the same way as
% the model or compatibility issues will arise (e.g., same preprocessing).
% Examples on how to use it:
%   Prediction = ClassificationPrediction(Data,Model_knn,'knn');
%   Prediction = ClassificationPrediction(Data,Model_LRegDA,'LReg-DA', ...
%       Options=Options_LRegDA);
%   Prediction = ClassificationPrediction(Data,Model_Kmeans,'Kmeans', ...
%       DoPCA=1,nComp=5);
% Please note that the syntax on how to specify option input has changed
% since Matlab 2021a. Example of before Matlab 2021a:
%   Prediction = ClassificationPrediction(Data,Model_Kmeans,'Kmeans',...
%       'DoPCA',1,'nComp',5);
% -------------------------------------------------------------------------
% Input:
%   Data:   The X-block of the data (the descriptors)
%               Its size will be m x n (m: number of samples; n: variables)
%   Model:  The model of the trained classification method
%   Method: The classification method
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
%   Options:    Options for 'PLS-DA' and 'LReg-DA' that can be provided.
%               See 'help lregda' or 'help plsda' for more information
%   DoPCA:      Do PCA before using KMeans or not (1/0 default: 1)
%   nComp:      The number of components to do the PCA with (only Kmeans).
%               Default: 5.
%
% Output:
%   Prediction: A column vector (m x 1) containing the class predictions of
%               the unknown samples.
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
DefaultOptions = [];
DefaultDoPCA = 1;
DefaultnComp = 5;

% Parse the input
p = inputParser;
addRequired(p,'Data',@ismatrix);
addRequired(p,'Model');
addRequired(p,'Method',@ischar);
addOptional(p,'Options',DefaultOptions);
addOptional(p,'DoPCA',DefaultDoPCA);
addOptional(p,'nComp',DefaultnComp);
parse(p,Data,Model,Method,varargin{:});

Data = p.Results.Data;
Model = p.Results.Model;
Method = p.Results.Method;
Options = p.Results.Options;
DoPCA = p.Results.DoPCA;
nComp = p.Results.nComp;

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

% Turn warnings off, to not flood the command window
warning('off','all');

% Check if the input Data is in the correct format.
if ~ismatrix(Data) || size(Data,2) == 1
    error("The input Data should be specified as a matrix (m x n: m = number of samples; n = number of variables).")
end

% Check if the specified method is part of the list
MethodCheck = ismember(Method, {'knn','RandomForest','PLS-DA','LReg-DA','AdaBoostM1','AdaBoostM2','RUSBoost','LogitBoost','HCA','Kmeans'});
if MethodCheck == 0
    error("The method can only be one of following (see 'help Classification' for more information): 'knn', 'RandomForest', 'PLS-DA', 'LReg-DA', 'AdaBoostM1', 'AdaBoostM2', 'RUSBoost', 'LogitBoost', 'HCA', 'Kmeans'")
end

% Predict the data using a trained classification model
switch Method

    % K-nearest neighbours classifier. This is a function that is built-in 
    % in the Statistics and Machine Learning toolbox of Matlab.
    case 'knn'

        % Autoscale the data, as it is not done within the knn routine
        Data = AutoScale(Data);

        % Predict the values of the data
        Prediction = predict(Model,Data);

    % Random Forest classifier. This is a function that is built-in 
    % in the Statistics and Machine Learning toolbox of Matlab.
    case 'RandomForest'

        % Autoscale the data, as it is not done within the Random Forest
        % routine
        Data = AutoScale(Data);

        % Predict the values of the data
        Prediction = predict(Model,Data);
        Prediction = cellfun(@str2num, Prediction);

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

        % Do the final model with an optimized number of components.
        ValidationModel = plsda(Data,Model, Options);
        Prediction = ValidationModel.classification.mostprobable;

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

        ValidationModel = lregda(Data,Model,Options);
        Prediction = ValidationModel.classification.mostprobable;

    % Adaptive boosting for Binary classifier. This is a function that 
    % is built-in in the Statistics and Machine Learning toolbox of Matlab.
    case 'AdaBoostM1'

        % Autoscale the data, as it is not done within the Random Forest
        % routine
        Data = AutoScale(Data);

        % Predict the values of the data
        Prediction = predict(Model,Data);

    % Adaptive boosting for multiclass classifier. This is a function that 
    % is built-in in the Statistics and Machine Learning toolbox of Matlab.
    case 'AdaBoostM2'

        % Autoscale the data, as it is not done within the Random Forest
        % routine
        Data = AutoScale(Data);

        % Predict the values of the data
        Prediction = predict(Model,Data);

    % Random Undersampling Boosting classifier. This is a function that 
    % is built-in in the Statistics and Machine Learning toolbox of Matlab.
    case 'RUSBoost'

        % Autoscale the data, as it is not done within the Random Forest
        % routine
        Data = AutoScale(Data);

        % Predict the values of the data
        Prediction = predict(Model,Data);

    % Adaptive Logistic Regression classifier. This is a function that is 
    % built-in in the Statistics and Machine Learning toolbox of Matlab.
    case 'LogitBoost'

        % Autoscale the data, as it is not done within the Random Forest
        % routine
        Data = AutoScale(Data);

        % Predict the values of the data
        Prediction = predict(Model,Data);

    % Hierarchical Cluster Analysis from the PLS toolbox. A 30-day trial 
    % demo can be downloaded from www.eigenvector.com
    case 'HCA'

        % No prediction can be done for hierarchical clustering on unknown
        % data.
        disp('No prediction can be done for the ''HCA'' method on unknown data. It is replaced by all 0s. A grouping is obtained from the training data.')
        Prediction = zeros(size(Data,1),1);

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

        % Predict the values of the data
        [~,Prediction] = pdist2(Model,Data,'euclidean','Smallest',1);

    % No variable selection. This should not occur as there is an error
    % check at the start of this function
    otherwise

        % Make an empty classification model
        Prediction = zeros(size(Data,1),1);

end

% Turn warnings back on
warning('on','all');

end