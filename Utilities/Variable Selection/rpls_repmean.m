function [Xmean,Xrep,Xstd]=rpls_repmean(X,rep)
% Finds the mean of the replicates (usually in predictions)
%[Xmean, Xrep, Xstd] = repmean(X,rep)
%
% INPUT:
%  X      X-matrix
%  rep    Sample index. Replicates have same number
%
% OUTPUT:
%  Xmean  Mean value of replicates
%  Xrep   Same length as X, but with average of the replicates
%  Xstd   Standard deviation of the replicates
%
% See also: MEAN

% Uses: fjernlike, length, ones, sparse, sum, zeros

% 150210 AAR Made the code work better with missing values
% 100207 AAR Added standard deviation
% 131106 AAR Adjusted it to also allow for those cases where there are no
%             replicates
% 141105 AAR Had made an error in the last correction
% 150805 AAR In some cases the rep-vector is not running. Make sure that
%             the algorithm also works in these cases.
% 230505 AAR Improved the algorithm considerably by using sparse matrix
%             notation

[i,j]=rpls_fjernlike(rep);
if length(i)==size(X,1)
    Xmean=X;
    Xrep=X;
else    
    %Here the replicates are recoded into increasing numbers
    rep = repnum( rep);%=vec((1:length(i))'*ones(1,size(j,2)));
%     j=vec(j);
%     rep(isnan(j))=[];
%     j(isnan(j))=[];
%     rep(j)=rep;
    
    s = sparse( 1:length( rep), rep, ones( length( rep), 1) );
    %Define a matrix which are all the NaN's of X
    nanX = isnan( X);
    %Calculate the sum of each set of replicates
    temp = X;
    temp( nanX ) = 0;
    temp = (temp' * s)';
    %Find the number of non NaN-numbers
    nanX = (~nanX' * s)';
    %Calculate the mean value
    Xmean = temp./ nanX;
    
    %Make the matrix which has the same size as X (but with mean values)
    Xrep = Xmean( rep, :);
    
    %Calculate the standard deviations
    temp = (X - Xrep).^2; %The deviations from the mean
    temp = (temp' * s)';
    nanX = nanX - 1;
    nanX( nanX == 0) = NaN;
    Xstd = sqrt( temp./ nanX);
end