function [Data,Class,NaNIdx,InfIdx] = RemoveNaNInf(Data,varargin)

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