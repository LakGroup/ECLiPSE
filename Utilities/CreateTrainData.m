function [TrainData,TrainClass,TestData,TestClass] = CreateTrainData(Data,Class,varargin)
% -------------------------------------------------------------------------
% Function that creates training data from a full data set. The full data
% set should contain multiple classes and then a random subset with the
% specified size (x number of classes) of the data will be selected as
% output.
% Examples on how to use it:
%   [TrainingData,Training_GT] = CreateTrainData(Data_PreProcessed,GroundTruth,ClassSize=[0.7 0.5]);
%   [TrainingData,Training_GT,ValidationData,ValidationClass] = CreateTrainData(Data_PreProcessed,GroundTruth,ClassSize=250,NeedTestData=1);
% -------------------------------------------------------------------------
% Input:
%   Data:           The X-block of the data (the descriptors)
%                       Its size will be m x n (m: number of samples; n:
%                       variables)
%   Class:          The Y-block of the data (the classes)
%                       Its size will be m x 1 (m: number of samples)
%
% Optional input:
%   ClassSize:      The number of samples per class in the training data
%   NeedTestData:   Whether or not test/validation data is needed.
%                   This is included in this function to make it
%                   independent of the training data
%
% Output:
%   TrainData:  The data that weas selected for the training dataset.
%   TrainClass: The classes of the data points inside the Traindata matrix.
%   TestData:   The data that was selected for the test/validation dataset.
%   TestClass:  The classes of the data points inside the TestData matrix.
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
DefaultClassSize = 0;
DefaultTest = 0;

% Parse the input
p = inputParser;
validScalar1 = @(x) isnumeric(x) && isscalar(x) && ~(x < 0) && ~(x > 1);
addRequired(p,'Data',@ismatrix);
addRequired(p,'Class',@ismatrix);
addOptional(p,'ClassSize',DefaultClassSize,@isnumeric);
addOptional(p,'NeedTestData',DefaultTest,validScalar1);
parse(p,Data,Class,varargin{:});

Data = p.Results.Data;
Class = p.Results.Class;
ClassSize = p.Results.ClassSize;
NeedTestData = p.Results.NeedTestData;

clear p DefaultTrainNumber DefaultTest validScalar1

% Check if the input Data is in the correct format.
if ~ismatrix(Data) || size(Data,2) == 1
    error('The input Data should be specified as a matrix (m x n: m = number of samples; n = number of variables).')
end

% Check if the input Class is in the correct format.
if ~ismatrix(Class) || size(Class,2) ~= 1
    error('The input Class should be specified as a vector (m x 1: m = number of samples).')
end

% Check the unique classes and check how many need to be sampled
[CountsClasses,UniqueClasses] = groupcounts(Class);
[minClass,Idx] = min(CountsClasses);

% Check if the input Class size is corretly specified.
if numel(ClassSize) > 1 && numel(ClassSize) ~= numel(UniqueClasses)
    error(['The input class size should be specified as a vector containing the same number of values as there are classes. Your data contains: ' numel(UniqueClasses) ' classes.'])
end
if numel(ClassSize) == 1 && ClassSize > min(CountsClasses)
    error(['The input class size should be lower than the number of samples in the smallest class (' num2str(minClass) ' for class ' num2str(UniqueClasses(Idx)) ')'])
end
if numel(ClassSize) > 1 && sum(gt(ClassSize,minClass)) >= 1
    error(['The input class size should be lower than the number of samples in each class: [' num2str(CountsClasses') ']'])
end

% Prepare the Train_Nmb & Test_Nmb vector for all classes
% If no size is specified, we go for a 70/30 split.
% If sizes are specified, we respect the sizes (fractions or whole numbers)
if numel(ClassSize) == 1
    if ClassSize == 0
        Train_Nmb = floor(minClass*0.7);
        Test_Nmb = min(CountsClasses) - Train_Nmb;
    else
        if ClassSize > 0 && ClassSize < 1
            Train_Nmb = floor(minClass*ClassSize);
            Test_Nmb = min(CountsClasses) - Train_Nmb;
        else
            Train_Nmb = ClassSize;
            Test_Nmb = minClass - Train_Nmb;
        end
    end
    Train_Nmb = repmat(Train_Nmb,[numel(UniqueClasses) 1]);
    Test_Nmb = repmat(Test_Nmb,[numel(UniqueClasses) 1]);
else
    if sum(ClassSize) == 0
        Train_Nmb = floor(CountsClasses*0.7);
        Test_Nmb = repmat(min(CountsClasses),[1 numel(CountsClasses)]) - Train_Nmb;
    else
        Train_Nmb = zeros(numel(UniqueClasses),1);
        Test_Nmb = zeros(numel(UniqueClasses),1);
        for i = 1:numel(ClassSize)
            if ClassSize(i) > 0 && ClassSize(i) < 1
                Train_Nmb(i) = floor(CountsClasses(i)*ClassSize(i));
                Test_Nmb(i) = min(CountsClasses) - Train_Nmb(i);
            else
                Train_Nmb(i) = ClassSize(i);
                Test_Nmb(i) = min(CountsClasses) - Train_Nmb(i);
            end
        end
    end
end

% Create the sampled data.
TrainData = zeros(sum(Train_Nmb),size(Data,2));
TrainClass = zeros(sum(Train_Nmb),1);

if NeedTestData == 1
    TestData = zeros(sum(Test_Nmb),size(Data,2));
    TestClass = zeros(sum(Test_Nmb),1);
end

% Loop over all the different classes and extract (without taking a same
% sample twice) the training data and training class
for i = 1:numel(UniqueClasses)
    ClassData = Data(Class==UniqueClasses(i),:);
    [TrainData(sum(Train_Nmb(1:i-1))+1:sum(Train_Nmb(1:i)),:),Idx] = datasample(ClassData,Train_Nmb(i),'Replace',false);
    TrainClass(sum(Train_Nmb(1:i-1))+1:sum(Train_Nmb(1:i)),:) = UniqueClasses(i);
    if NeedTestData == 1
        Idx = setdiff(1:size(ClassData,1),Idx);
        ClassDataTest = ClassData(Idx,:);
        TestData(sum(Test_Nmb(1:i-1))+1:sum(Test_Nmb(1:i)),:) = datasample(ClassDataTest,Test_Nmb(i),'Replace',false);
        TestClass(sum(Test_Nmb(1:i-1))+1:sum(Test_Nmb(1:i))) = UniqueClasses(i);
    end
end

if NeedTestData ~=1
    TestData = [];
    TestClass = [];
end

end