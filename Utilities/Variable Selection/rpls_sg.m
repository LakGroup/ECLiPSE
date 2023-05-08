function xnew=rpls_sg( x, m, d, k, opt, incl)
%Smooths and/ or calculates the derivative according to Savitsky-Golay
%
%xnew=rpls_sg(x,m,d,k,opt, incl)
%
%INPUT:
% x    The raw spectra. Samples in rows
% m    Number of points to be used for the smoothing (2m+1, m=1:15) <Default = 3>
% d    Degreee of derivative. d=0:5 <Default = 1>
% k    Degree of polynom to fit to the curve. k=2:6 <Default = 3> 
% opt  Inclusion of end-points (1) or not (0) <Default = 1>
% incl What points to include in the smoothing. 
%       <Default = true( 2m + 1, 1)>
%
%OUTPUT:
% xnew The smoothed (and/or derivative)
%
% See also: MSC, SNV, SPECDET, SGCONV

% 211216 AAR Included the possibility to not include all points (i.e. for
%             the removal of spikes
% 301007 AAR Fixed the calculations of the end-points
% 081007 AAR Added the end-points
% 250907 AAR Increased the possible values for m, d and k
% 090107 AAR

if nargin==0
    error('No input arguments given')
end
%Set default values
if nargin<2 || m<1 %|| m>15
    m=3;
    disp('''m'' set to 3')    
end
if nargin<3 || d<0 || d>4
    d=1;
    disp('''d'' set to 1')
end
if nargin<4 || k<0 || k>8
    k=3;
    disp('''k'' set to 3')
end
if nargin<5 || opt>1 || opt<0
    opt=1;
end

if nargin < 6
    incl = true( 2 * m + 1, 1);
end

if size(x,2)==1
    error('The samples should be in the rows, not the columns')
end

convolute=rpls_sgconv( m, d, k, 0, incl);

if isempty(convolute)
    error(['The combination m=' num2str(m) ', d=' num2str(d) ' and k=' num2str(k) ' is not possible'])
end
xnew=convn(x,convolute,'valid');

xnew=[ones(size(xnew,1),m)*NaN xnew ones(size(xnew,1),m)*NaN];

if opt==1
    %Inclusion of the end-points
    %Reference Gorry, Analytical Chemistry, 1990, 62, 570-573.
    for i=1:m
        convs=rpls_sgconv(m,d,k,i, incl);
        conve=rpls_sgconv(m,d,k,-i, incl);
        xnew(:,(m+1)-i)=convn(x(:,1:2*m+1),convs,'valid');
        xnew(:,end-(m-i))=convn(x(:,end-2*m:end),conve,'valid');
    end
end