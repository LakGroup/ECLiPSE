function a=rpls_nansum(x, dim)
% a=rpls_nansum(x, dim)
%
% Calculate the sum without taking into account the missing values.
%
% INPUT
%  x    X-matrix/vector
%  dim  The dimension along which the sum is calculated < default = 1 >
%
% OUTPUT
%  a    Sum-values
%
% See also: nanmax, nanmean, nanmedian, nanmin, nanstd, nansum, nanvar

% 221008 AAR Allowed for calculating the sum in any direction
% 071207 AAR Corrected for situations where a column has all NaN's
% 010905 AAR

if nargin < 2
    dim = 1;
end
if dim > length( size( x) )
    dim = 1;
end

c = x;
x(isnan(x))=0;
a = sum(x, dim);
siz = size(x);
a( sum( isnan(c), dim) == siz(dim)) = NaN;