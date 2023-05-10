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
%   Hugelier, S., Kim, H., Gyparaki, M.T., Bond, C., Tang, Q., 
%   Santiago-Ruiz, A.N., Porta, S., Lakadamyali, M. ECLiPSE: a versatile 
%   classification technique for structural and morphological analysis of 
%   super-resolution microscopy data. BioRxiv (2023). 
%   DOI: https://doi.org/10.1101/2023.05.10.540077.
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