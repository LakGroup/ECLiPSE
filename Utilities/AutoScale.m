function AutoScaledData = AutoScale(Data)
% -------------------------------------------------------------------------
% Function that autoscales the data. Autoscalimg is unit variance scaling.
% Example on how to use it:
%   Data_Preprocessed = AutoScale(RawData);
% -------------------------------------------------------------------------
% Input:
%   Data:   The data to be autoscaled
%               Its size should be m x n (m: number of samples; n: variables)
%
% Output:
%   AutoScaledData: The autoscaled data
% -------------------------------------------------------------------------
% Code written by:
%   Siewert Hugelier    Lakadamyali lab, University of Pennsylvania (USA)
% Contact:
%   siewert.hugelier@pennmedicine.upenn.edu
%   melike.lakadamyali@pennmedicine.upenn.edu
% If used, please cite:
%   xxx
% -------------------------------------------------------------------------

% Check if the input Data is in the correct format.
if ~ismatrix(Data) || size(Data,2) == 1
    error("The input Data should be specified as a matrix (m x n: m = number of samples; n = number of variables).")
end

% Calculate the mean and standard deviation of the data
meanData = mean(Data);
stdData = std(Data);

% Calculate the autoscaled data. Autoscaling data is also known as unit
% variance scaling. It consists of two steps: (1) mean centering, and (2)
% scaling/standardization
AutoScaledData = (Data - meanData(ones(size(Data,1),1),:)) ./ stdData(ones(size(Data,1),1),:);

end