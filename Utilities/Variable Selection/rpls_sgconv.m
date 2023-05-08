function conv=rpls_sgconv( m, d, k, i, incl)
%Calculate the convolutes for Savitsky-Golay derivation (d) for 2m+1
%smoothing with a k-degree polynomial for point i (out of 2m+1)
% conv = rpls_sgconv( m, d, k, i, incl)
%
%INPUT:
% m     2m+1 points for smoothing
% d     Derivative (d>=0)
% k     Degree of polynomial fitting across the 2m+1 points
% i     Convolution point (Default = 0)
% incl  What points should be included in the fitting
%        <Default = true( 2m + 1, 1) >
%
%OUTPUT:
% conv  The convolute

% 081007 AAR Allowed for all i's
% 250907 AAR

%Set the default value
if nargin<4 || i>m || i<-m
    i=0;
end

if nargin < 5
    incl = true( 2 * m + 1, 1);
end

c=zeros(2*m+1,k);
c(:,1)=((m:-1:-m)')+i;
for i=1:k
    c(:,i)=c(:,1).^i;
end
c=[ones(length(c),1) c];
c = c( incl, :);
if rank(c) < size(c,2)
    temp = pinv( c, eps*3);
else
    temp = inv(c'*c)*c';
end

conv = zeros( 1, length( incl) );
%In case something illegal is asked for

if (d+1) > size(temp,1)
    conv = [];
else
    conv( incl) = temp( d+1, :) * factorial(d);
end
