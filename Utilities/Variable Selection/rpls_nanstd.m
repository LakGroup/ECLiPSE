function a=rpls_nanstd(x, dim)
% a=rpls_nanstd(x, dim)
%
% Calculate the var without taking into account the missing values.
%
% INPUT
%  x    X-matrix/vector
%  dim  The dimension the std should be caclulated for (default = 1)
%
% OUTPUT
%  a    Std-values
%
% See also: nanmax, nanmean, nanmedian, nanmin, nanstd, nansum, nanvar

% 090211 AAR Make it work with dimensions
% 010905 AAR

if nargin == 1
    dim = 1;
end

if dim > 2 || dim < 1
    dim = 1;
    disp( '''dim'' has been set to 1')
end

m=rpls_nanmean(x);
x=x-ones(size(x,1),1)*m;
x=x.^2;
a = sqrt( rpls_nansum( x, dim)./ (sum( ~isnan( x), dim) - 1) );