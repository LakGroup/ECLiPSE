function [Data,Class,NaNIdx,InfIdx] = RemoveNaNInf(Data,varargin)
% -------------------------------------------------------------------------
% Function that removes any rows or columns from the data that contain NaNs
% or Inf values.
% Examples on how to use it:
%   [Data,Class,NaNIdx,InfIdx] = RemoveNaNInf(Data,Class=GT,Col=0);
%   Data = RemoveNaNInf(Data,Col=0);
% Please note that the syntax on how to specify option input has changed
% since Matlab 2021a. Example of before Matlab 2021a:
%   [Data,Class,NaNIdx,InfIdx] = RemoveNaNInf(Data,'Class',GT,'Col',0);
% -------------------------------------------------------------------------
% Input:
%   Data:   The data matrix of the descriptors in which the Inf and NaN
%           values have to be removed.
%           Its size is m x n (m: number of samples, n: number of
%           variables).
%
% Optional input:
%   Class:  A vector with known class associations of the samples in the
%           data matrix.
%           Its size is m x 1 (m: number of samples).
%   Col:    A scalar representing whether rows (0) or columns (1) have to
%           be removed. Default: 0.
%
% Output:
%   Data:   The data matrix in which the rows/columns containing NaN and/or
%           Inf values have been removed.
%   Class:  The vector containing the known class associations in which the
%           entries with NaN or Inf values are removed.
%   NaNIdx: The indices of the rows/colums containing in NaN values.
%   InfIdx: The indices of the rows/colums containing in Inf values.
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
DefaultClass = [];
DefaultCol = 0;

% Parse the input
p = inputParser;
addRequired(p,'Data',@ismatrix);
addOptional(p,'Class',DefaultClass);
addOptional(p,'Col',DefaultCol);
parse(p,Data,varargin{:});

Data = p.Results.Data;
Class = p.Results.Class;
Col = p.Results.Col;

% Check the input data and convert to a matrix if it is a table
if istable(Data)
    Data = table2array(Data);
end

% Do it either on rows or columns
if Col ~= 1
    % Check for NaN values in the data
    NaNIdx = any(isnan(Data'));
    Data(NaNIdx,:) = [];
    
    if ~isempty(Class)
        Class(NaNIdx) = [];
    end
    
    % Check for Inf values in the data
    InfIdx = any(isinf(Data),2);
    Data(InfIdx,:) = [];
    
    if ~isempty(Class)
        Class(InfIdx,:) = [];
    end
    
else
    % Check for NaN values in the data
    NaNIdx = any(isnan(Data));
    Data(:,NaNIdx) = [];
    
    % Check for Inf values in the data
    InfIdx = any(isinf(Data),1);
    Data(:,InfIdx) = [];
end

if nargout >= 2 && isempty(Class)
    Class = ones(1,size(Data,1));
end

end