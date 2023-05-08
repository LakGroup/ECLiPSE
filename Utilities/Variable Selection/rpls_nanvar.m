function xvar=rpls_nanvar(x, dim)
%xvar=rpls_nanvar(x, dim)
%
% Variance ignoring NaNs
%
% INPUT:
%  x     X-matrix
%  dim   The dimension the std should be caclulated for (default = 1)
%
% OUTPUT:
%  xvar  Variance of the columns in x, ignoring the NaNs
%
% See also: nanmax, nanmean, nanmedian, nanstd, nansum

%AAR 310804

if nargin < 2
    dim = 1;
end

if dim > length( size( x) )
    error( 'Your data does not have such a high dimension')
end

xvar = rpls_nanstd( x, dim).^2;