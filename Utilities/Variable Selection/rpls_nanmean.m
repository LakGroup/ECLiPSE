function a=rpls_nanmean(x, dim)
% a=rpls_nanmean(x, dim)
%
% Calculate the mean without taking into account the missing values.
%
% INPUT
%  x    X-matrix/vector
%  dim  Direction. Default = 1 (column-wise), 2 = row-wise
%
% OUTPUT
%  a    Mean-values
%
% See also: NANMAX, NANMEDIAN, NANMIN, NANSTD, NANSUM, NANVAR

% 090211 AAR Make the output similar to mean for dimensions
% 250111 AAR Made it work with dimensions
% 241105 AAR 
% 010905 AAR

x = x;
if nargin == 1
    if min( size( x) ) == 1
        x = rpls_vec( x);
    end
    dim = 1;
end

if dim < 1 || dim > 2
    error( '''dim'' should be either ''1'' or ''2''')
end

c=~isnan(x);
x(isnan(x))=0;

sumx=sum(x, dim);
sumc=sum(c, dim);
a(sumc>0)=sumx(sumc>0)./sumc(sumc>0);
a(sumc==0)=NaN;
if dim == 2
    a = a';
end