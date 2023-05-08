function VarSelModel = VariableSelection(Data,Class,Method)
% -------------------------------------------------------------------------
% Function that uses automatic variable/feature selection to focus on the
% important variables for e.g., classification. Depending on the method
% used, the results will select more or less variables. A selection of
% different methods can also be used and then combined to select variables
% that appear multiple times.
% Note that some toolboxes need to be used for some of the methods to work
% (see below). The iToolbox and the rPLS toolbox are included in the
% functions. The PLS Toolbox and the Statistics & Machine Learning Toolbox
% are not included.
% Examples on how to use it:
%   VarSel_iPLS = VariableSelection(TrainingData1,Training_GT1,'iPLS_PLSToolbox');
%   VarSel_biPLS = VariableSelection(TrainingData2,Training_GT2,'biPLS');
% -------------------------------------------------------------------------
% Input:
%   Data:   The X-block of the data (the descriptors)
%               Its size will be m x n (m: number of samples; n: variables)
%   Class:  The Y-block of the data (the classes)
%               Its size will be m x 1 (m: number of samples)
%   Method: The variable selection method
%           Choices:    'biPLS' (Backwards interval PLS; iToolbox)
%                       'rPLS' (Recursive PLS; kindly provided by Prof. Åsmund Rinnen)
%                       'biPLS_PLSToolbox' (Interval PLS from the PLS Toolbox; Eigenvector Inc.)
%                       'rPLS_PLSToolbox' (Recursive PLS from the PLS Toolbox; Eigenvector Inc.)
%                       'iPLS_PLSToolbox' (Interval PLS from the PLS Toolbox; Eigenvector Inc.)
%                       'GA_PLSToolbox' (Genetic Algorithm from the PLS Toolbox; Eigenvector Inc.)
%                       'ChiSquare' (Chi Square test from the Statistics & Machine Learning Toolbox; The MathWorks Inc.)
%                       'MRMR' (Minimum Redundancy Maximum Relevance algorithm from the Statistics & Machine Learning Toolbox; The MathWorks Inc.)
%                       'ReliefF' (RliefF algorithm from the Statistics & Machine Learning Toolbox; The MathWorks Inc.)
%                       'Boruta' (Boruta algorithm - simplistic implementation)
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

% Update the path to make sure the PLS toolbox is in the appropriate
% position, depending on the method being used
if contains(Method,'PLSToolbox')
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

% Check if the input Class is in the correct format.
if ~ismatrix(Class) || size(Class,2) ~= 1
    error("The input Class should be specified as a vector (m x 1: m = number of samples).")
end

% Check if the specified method is part of the list
MethodCheck = ismember(Method, {'biPLS','rPLS','biPLS_PLSToolbox','rPLS_PLSToolbox','iPLS_PLSToolbox','GA_PLSToolbox','ChiSquare','MRMR','ReliefF','Boruta'});
if MethodCheck == 0
    error("The method can only be one of following (see 'help VariableSelection' for more information): 'biPLS','rPLS','biPLS_PLSToolbox','rPLS_PLSToolbox','iPLS_PLSToolbox','GA_PLSToolbox','ChiSquare','MRMR','ReliefF','Boruta'")
end

% Perform the variable selection
switch Method

    % Backwards interval partial least squares: slower but more
    % accurate according to our experience. Comes from the iToolbox:
    % http://www.models.life.ku.dk/itoolbox
    case 'biPLS'

        % Initialization
        VarSel = zeros(1,size(Data,2));

        % Set up the model
        model = bipls_NotPLSToolbox(Data,Class,20,'auto',size(Data,2),[],'syst123',5);

        % Extract the variables that were kept for optimal models
        Intervals = zeros(size(model.allint,1),max(model.allint(1:end-1,3)+1-model.allint(1:end-1,2)));
        for j = 1:size(model.allint,1)-1
            Intervals(j,1:model.allint(j,3)+1-model.allint(j,2)) = (model.allint(j,2):model.allint(j,3));
        end

        [~,Idx] = min(model.RevRMSE); % Check which model was the most accurate

        Intervals(model.RevIntInfo(1:Idx-1),:) = [];
        Intervals = Intervals(:);
        Intervals(Intervals==0) = [];
        VarSel(Intervals) = 1;

        % Create output
        VarSelModel.model = model;
        VarSelModel.VarSel = VarSel;
        VarSelModel.Method = Method;

    % Recursive partial least squares: faster but gives less
    % accurate results in the applications we tested. It seems like
    % this is more useful when variables are continuous
    % (spectroscopy). Kindly provided by Prof. Åsmund Rinnen
    case 'rPLS'

        % Autoscale the data, as it is not done within the rPLS routine
        Data = AutoScale(Data);

        % Set all the options needed for the cross validation (to
        % calculate optimal models)
        options = rpls_NotPLSToolbox;
        options.nLV = 8;
        options.Set = [Inf 1e14];
        options.Val.ind = rpls_cvind(Data, [], 10);

        % Perform the variable selection. We do a check here
        % because sometimes it gives an unexplicable error, and so
        % we just ignore it and continue with the variable selection.
        try
            [model, evolve] = rpls_NotPLSToolbox(Data, Class, options);

            % Extract the variables that were kept for optimal models
            VarSel = double(evolve.ind(evolve.optimal,:));
        catch
            model = [];
            evolve = [];
            VarSel = [];
        end

        % Create output
        VarSelModel.model = model;
        VarSelModel.evolve = evolve; % The only one that needs this
        VarSelModel.VarSel = VarSel;
        VarSelModel.Method = Method;

    % Backwards iPLS from the PLS toolbox. A 30-day trial demo can
    % be downloaded from www.eigenvector.com
    case 'biPLS_PLSToolbox'

        % Initialization
        VarSel = zeros(1,size(Data,2));

        % Create the options
        options = ipls('options');
        options.preprocessing = {'autoscale' 'autoscale'};
        options.mode = 'reverse';
        options.plots = 'none';
        options.display = 'off';

        % Perform the variable selection
        model = ipls(Data,Class,1,round(size(Data,2)/5),Inf,options);

        % Extract the variables that were kept for optimal model
        VarSel(model.use) = 1;

        % Create output
        VarSelModel.model = model;
        VarSelModel.VarSel = VarSel;
        VarSelModel.Method = Method;

    % Recursive iPLS from the PLS toolbox. A 30-day trial demo can
    % be downloaded from www.eigenvector.com
    case 'rPLS_PLSToolbox'

        % Initialization
        VarSel = zeros(1,size(Data,2));

        % Create the options
        options = rpls('options');
        options.mode = 'suggested';
        options.preprocessing = 'autoscale';
        options.plots = 'none';
        options.display = 'off';
        options.waitbar = 'off';

        % Perform the variable selection
        model = rpls(Data,Class,round(size(Data,2)/5),options);

        % Extract the variables that were kept for optimal model
        [~,Idx] = min(model.RMSECVs);
        VarSel(model.selectedIdxs{Idx}) = 1;

        % Create output
        VarSelModel.model = model;
        VarSelModel.VarSel = VarSel;
        VarSelModel.Method = Method;

    % iPLS from the PLS toolbox. A 30-day trial demo can be
    % downloaded from www.eigenvector.com
    case 'iPLS_PLSToolbox'

        % Initialization
        VarSel = zeros(1,size(Data,2));

        % Create the options
        options = ipls('options');
        options.preprocessing = {'autoscale' 'autoscale'};
        options.plots = 'none';
        options.display = 'off';

        % Perform the variable selection
        model = ipls(Data,Class,1,round(size(Data,2)/5),Inf,options);

        % Extract the variables that were kept for the model
        VarSel(model.use) = 1;

        % Create output
        VarSelModel.model = model;
        VarSelModel.VarSel = VarSel;
        VarSelModel.Method = Method;

    % Genetic Algorithm from the PLS toolbox. A 30-day trial demo
    % can be downloaded from www.eigenvector.com
    case 'GA_PLSToolbox'

        % Create the options
        options = gaselctr('options');
        options.preprocessing = {preprocess('default','autoscale') []};
        options.plots = 'none';
        options.display = 'off';
        options.popsize = 64; % More takes longer, but more representative of the data
        options.reps = 10; % Best to use a low amount of generations and more repetitions (avoids overfitting)
        options.maxgenerations = 10; % Best to use a low amount of generations and more repetitions (avoids overfitting)
        options.cv = 'rnd'; % Best for genetic algorithms

        % Perform the variable selection
        model = gaselctr(Data,Class,options);

        % Extract the variables that were kept for optimal model
        [~,Idx] = min(model.rmsecv);
        VarSel = model.icol(Idx,:);

        % Create output
        VarSelModel.model = model;
        VarSelModel.VarSel = VarSel;
        VarSelModel.Method = Method;

        % Feature selection based on Chi Square tests. This needs the
        % statistics and machine learning toolbox of Matlab
    case 'ChiSquare'

        % Autoscale the data, as it is not done within the fscchi2 routine
        Data = AutoScale(Data);

        % Perform the variable selection
        [Idx,Scores] = fscchi2(Data,Class);

        % Extract the variables that were important for the model
        [~,Idx] = sort(Idx);
        if any(isinf(Scores))
            VarSel = zeros(1,size(Data,2));
            Scores = Scores(Idx);
            VarSel(isinf(Scores)) = 1;
        else
            VarSel = double(Scores(Idx)./sum(Scores)*100 > 100/size(Data,2)); % More important than 1 single variable
        end

        % Create output
        VarSelModel.model = [];
        VarSelModel.VarSel = VarSel;
        VarSelModel.Method = Method;

    % Feature selection based on the Minimum Redundancy Maximum
    % Relevance (MRMR) Algorithm. This needs the statistics and
    % machine learning toolbox of Matlab
    case 'MRMR'

        % Autoscale the data, as it is not done within the fscmrmr routine
        Data = AutoScale(Data);

        % Perform the variable selection
        [Idx,Scores] = fscmrmr(array2table(Data),Class);

        % Extract the variables that were important for the model
        [~,Idx] = sort(Idx);
        VarSel = double(Scores(Idx) > 0.1);

        % Create output
        VarSelModel.model = [];
        VarSelModel.VarSel = VarSel;
        VarSelModel.Method = Method;

    % Feature selection based on the Minimum Redundancy Maximum
    % Relevance (MRMR) Algorithm. This needs the statistics and
    % machine learning toolbox of Matlab
    case 'ReliefF'

        % Initialize
        numNeighbours = 30; % 30 comes from a study that was done before for our application

        % Autoscale the data, as it is not done within the fscmrmr routine
        Data = AutoScale(Data);

        % Perform the variable selection
        [Idx,Scores] = relieff(Data,Class,numNeighbours,'method','classification');

        % Extract the variables that were kept for optimal model
        [~,Idx] = sort(Idx);
        VarSel = double(Scores(Idx)./sum(Scores)*100 > 100/size(Data,2)); % More important than 1 single variable

        % Create output
        VarSelModel.model = [];
        VarSelModel.VarSel = VarSel;
        VarSelModel.Method = Method;

    % Simplistic implementation of the Boruta R package, which can
    % be found at: 10.18637/jss.v036.i11
    case 'Boruta'

        % Initialization
        Confidence = 0.999; % This can be changed. Just standard value from the Boruta paper
        BorutaRuns = ceil(-log2(1-Confidence));
        HigherPct = binoinv(Confidence,BorutaRuns,0.5);
        LowerPct = binoinv(1-Confidence,BorutaRuns,0.5);

        % Start the waitbar
        wb = waitbar(0,['Performing Boruta Variable Selection. Please be patient... At iteration: 0/' num2str(BorutaRuns)]);

        % Autoscale the data, as it is not done within the routine
        Data = AutoScale(Data);

        BorutaCount = zeros(BorutaRuns,size(Data,2));
        % Perform the Boruta code
        for j = 1:BorutaRuns

            % Update the waitbar
            waitbar(j/BorutaRuns,wb,['Performing Boruta Variable Selection. Please be patient... At iteration: ' num2str(j) '/' num2str(BorutaRuns)]);

            % Initialize the Boruta training data.
            BorutaData = [Data zeros(size(Data,1),size(Data,2))];

            % Make the training data. Randomize each column of the
            % training data and add it as the second half of the
            % data
            for k = 1:size(Data,2)/2
                RandIdx = randperm(size(Data,1));
                BorutaData(:,size(Data,2)+k) = Data(RandIdx,k);
            end

            % Do the Boruta variable selection
            Model_Boruta = TreeBagger(100,BorutaData,Class,'Method','classification','OOBPredictorImportance','On','PredictorSelection','Curvature');
            BorutaImportance = Model_Boruta.OOBPermutedPredictorDeltaError;

            MaxImportance = max(BorutaImportance(size(Data,2)+1:end));
            BorutaCount(j,:) = double((BorutaImportance(1:size(Data,2)) > MaxImportance));

            clear k RandIdx BorutaData
        end
        close(wb)

        % Check which variables are definitely important, and
        % definitely not important
        model.BorutaCount = sum(BorutaCount);
        model.Support = double(model.BorutaCount > HigherPct);
        model.WeakSupport = double(model.BorutaCount >= LowerPct & model.BorutaCount <= HigherPct);
        model.NoSupport = double(model.BorutaCount < LowerPct);

        % Create output
        VarSelModel.model = model;
        VarSelModel.VarSel = model.Support;
        VarSelModel.Method = Method;

    % No variable selection. This should not occur as there is an error
    % check at the start of this function
    otherwise

        % Set all variables to 'important'
        VarSel = ones(1,size(Data,2));

        % Create output
        VarSelModel.model = [];
        VarSelModel.VarSel = VarSel;
        VarSelModel.Method = Method;

end

% Turn warnings back on
warning('on','all');

end