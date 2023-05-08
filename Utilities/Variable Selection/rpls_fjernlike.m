function [svar,ind] = rpls_fjernlike(x,n,lim)
%[svar,ind]=fjernlike(x,n,lim)
%
% Removes equal rows with respect to a specific column. NaN-values will not
% be counted as a unique value
%
% INPUT:
%  x      The x-matrix, can also be a string matrix
%  n      The column to check by 
%  lim    Can be added to specify a %-deviation allowed for equality
%
% OUTPUT:
%  svar   The matrix consisting of the unique values in column n
%  ind    The indexes for the equal rows

% Copyright, 2005 - 
% This M-file and the code in it belongs to the holder of the
% copyrights and is made public under the following constraints:
% It must not be changed or modified and code cannot be added.
% The file must be regarded as read-only. 
%
% Åsmund Rinnan
% Quality and Technology
% Department of Food Science
% Faculty of Life Sciences
% University of Copenhagen
% Rolighedsvej 30, 1958 Frederiksberg C, Denmark
% Phone:  +45 35 33 35 42
% e-mail: aar@life.ku.dk

% AAR 250413 Now also allows for cell arrays with several columns. The
%             output is the same format as the input (i.e. input cells give output
%             cells)
% AAR 140507 Now also handles character vectors
% AAR 271006 Gives now the correct ind, so that x(ind(12,:),:)=svar(12,:)
% AAR 300605 Made the algorithm notably faster for both numbers and
%             strings.
% AAR 290605 Made the algorithm faster by checking for the length of the
%             unique and similar components. (Can mean a lot if the number
%             of rows are high!) (Not implemented for strings)
% AAR 290605 Change cell to char
% AAR 170105 Changed if 'x' is string to use b~=
% AAR 130105 Changed 'ind' so that it is not made unnecessarily large

if isempty(x)
    svar=[];
    ind=[];
    return
end
if iscell(x)
    [rx, cx] = size( x);
    if size( x, 2) > 1
        test = 2;
        txt = char( x( :, 1) );
        txt( :, end + 1:end + 2) = ones( rx, 1) * '; ';
        for cs = 2:size( x, 2)
            txt = [txt char( x( :, cs) )];
            txt( :, end + 1:end + 2) = ones( rx, 1) * '; ';
        end
        x = txt( :, 1:end-2);
    else
        test = 1;
        x=char(x);
    end
else
    test = 0;
end

if nargin<3
    lim=0;
end
[r,k]=size(x);
if nargin==1
    n=1:k;
end
if isstr(x(1,:))
    [xnew,i]=sortrows(x,n);
    a=double(xnew);
    c=a(1:end-1,:)~=a(2:end,:);
    d=a(1:end-1,:)==a(2:end,:);
    %If there is only a character vector, this will also check for those
    if size(c,2)>1
        b=[1 find(sum(c')>0)+1]';
    else
        b=[1;find(c>0)+1]; 
    end
else
    [xnew,i]=sortrows(x(:,n));
    a=xnew(1:end-1,:)*(1+lim)-xnew(2:end,:)*(1-lim);
    a(a>0)=0;
    a=a.^2;
    if length(n)>1
        a=sum(a')';
    end
    b=[1;find(a~=0)+1];
end
svar=x(sort(i(b)),:);
if ~isstr(x(1,:))
    if size(svar,2)==1
        svar=svar(~isnan(svar));
    else
        svar=svar(sum(~isnan(svar)')==size(svar,2),:);
    end
end

if nargout==2
    %Find the unique rows
    a=zeros(size(x,1),1);
    a(b)=1;
    c=cumsum(a);
    %Find the equal rows
    a(a==1)=b-[0;b(1:end-1)];
    d=(1:length(a))'-cumsum(a)+1;
    [o,p]=sort(i(b));
    %Convert that into a row index with original coordinates
    ind=full(spconvert([c d i]));
    ind(ind==0)=NaN;
    %It is utterly important to include the following line!
    %If not the index of row 12 will NOT(!) be the unique row given in svar
    ind=ind(p,:);
end

switch test
    case 1 %Single column cells
        svar = cellstr( svar);
    case 2 %Multi-column cells
        id = [-1 strfind( svar( 1, :), ';') size( svar, 2) + 1];
        clear temp
        for cc = 1:cx
            temp( :, cc) = cellstr( svar( :, id(cc) + 2:id( cc + 1) - 1) );
        end
        svar = temp;
end