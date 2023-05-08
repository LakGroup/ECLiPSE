function [b,w,t,p, u, q] = rpls_bipls( x, y, nLV, nrm)
%[b, w, t, p, u, q] = bipls( x, y, nLV, nrm)
% Fitting the PLS by using the bidiagonal algorithm
%
%INPUT:
% x    X-matrix
% y    y-vector
% nLV  Number of factors to estimate
% nmr  Normalize the loadings to length 1 (default = false)
%
%OUTPUT:
% b    Regression coefficients
% w    Loading weights
% t    X-scores
% p    X-loadings
% u    y-scores
% q    y-loadings

%Based on: Martin Andersson: A comparison of nine PLS1 algorithms, Journal
%           of Chemometrics 2009; 23: 518-529
%Adjusted according to Ulf Indahl (unpublished work)
 
% Copyright, 2014 - 
% This M-file and the code in it belongs to the holder of the
% copyrights and is made public under the following constraints:
% It must not be changed or modified and code cannot be added.
% The file must be regarded as read-only. 
% In case of doubt, contact the holder of the copyrights.
%
% Åsmund Rinnan
% E-mail asmundrinnan@gmail.com

if nargin < 4
    nrm = false;
end

w(:,1) = x'*y;
w(:,1) = w(:,1)/ sqrt( w(:,1)' * w( :, 1) );
t(:,1) = x*w(:,1);
B(1,1) = sqrt( t( :, 1)' * t( :, 1) );
t(:,1) = t(:,1)/B(1,1);
nLV = min( rank( x), nLV);
for a = 2:nLV
    w(:,a) = x'*t(:,a-1)-B(a-1,a-1)*w(:,a-1);
    
    %Need to re-orthogonalize w
    corr =  w( :, a)' * w( :, 1:a-1);
    w( :, a) = w( :, a) - w( :, 1:a-1) * corr'; 

    B(a-1,a) = sqrt( w( :, a)' * w( :, a) );
    w(:,a) = w(:,a)/B(a-1,a);
    t(:,a) = x*w(:,a)-B(a-1,a)*t(:,a-1);
    
    %Need to re-orthogonalize t
    corr = t( :, a)' * t( :, 1:a-1);
    t( :, a) = t( :, a) - t( :, 1:a-1) * corr';
    
    B(a,a) = sqrt( t( :, a)' * t( :, a) );
    t(:,a) = t(:,a)/B(a,a);
end
invB = pinv(B);
q = y'*t;

% %A rough reduction of the nLVs
% bd = diag( B);
% temp = bd/ bd(1);
% nLV = min( find( temp < sqrt( eps), 1), nLV);

%Put the length of t back into t (as is normally done), and calculate p
told = t;
t = t * diag( diag( B));
p = ( pinv( t) * x)';

for a = 1:nLV
    b(:,a) = w(:,1:a)*(invB(1:a,1:a)*q(1:a)');
    u( :, a) = y/ q(a);
    y = y - told( :, a) * q(a);
end

if nrm
    t=t*norm(p);
    w=w*norm(p);
    p=p/norm(p);
end
