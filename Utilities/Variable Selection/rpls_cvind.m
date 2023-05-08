function ind = rpls_cvind(Y,rep,seg, ran)
%ind = rpls_cvind( Y, rep, seg, ran)
%
% INPUT:
%  Y      Values of the dependent variable
%  rep    Vector (numbers or string) representing the samplenames
%          (/-numbers) Default = 1:length(y)
%  seg    Number of segments for the cross-validation, or %/part of set in each 
%          segment (normalized to sum 1) -> for making Cal+(Stop)+Test
%  ran    Two number:
%          1 - Partly randomization of the index (True), still making sure that the
%               response space is covered nicely. Default = false
%          2 - The numbers are NOT classes. Default = true (i.e. NOT class)
%
% OUTPUT:
%  ind    Index for the crossvalidation segments
%
% See also: cv

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
    ran = false;
end

if isempty(rep)
    rep = rpls_vec(1:length(Y));
end
rep = rpls_vec(rep);
[ m, n] = rpls_fjernlike( rep);
[ i, j] = sort( Y( n( :, 1) ) ); %Only include one from each replicate
rep = zeros( length( rep), 1);
rep( n( :, 1) )= 1:length( n( :, 1) );
while any( rep == 0)
    zero = find( rep==0);
    rep( zero) = rep( zero-1);
end
r=max(rep);
if length(seg)>1
    per=seg./sum(seg);
    seg=length(i);
else
    per=0;
end
ind=zeros(size(rep));

if per == 0
    
    %150610 AAR Added lines for randomization of segmentation
    if ran(1)
        cvpos = zeros( length( m), 1);
        %A rough check to see if there are classes
        if length( ran) == 2
            test = ran( 2);
        else
            if length( m) < length( rpls_fjernlike( i))*3 %Increased this to 3
                test = true;
            else
                test = false;
            end
        end
        if test
            %Grouping the similar Y-values
            num = floor( length( m)/ seg);
            temp = ones( ceil( length( m)/num), 1) * (1: num);
            tg = prod( size( temp) ) - length( m);
            if tg > 0
                [tf, tl] = sort( rand( size( temp, 2), 1) );
                temp( end, tl( 1:tg) ) = NaN;
                temp = rpls_vec( temp);
                temp( isnan( temp) ) = [];
            else
                temp = rpls_vec( temp);
            end
            class( j, 1) = temp;
            [cl, cb] = rpls_fjernlike( class);
        else
            [cl, cb] = rpls_fjernlike( i); %Find the number of samples in each class
        end
        
        %To get the total number in each segment (from each group)
        ng = sum( ~isnan( cb') );
        ng = floor( ng/ seg);
        for cs = 1:length( cl)
            %I distribute the samples which should be in different
            %segments according to 'ng'
            [temp, id] = sort( rand( sum( ~isnan( cb( cs, :) ) ), 1) );
            id( ( ng(cs) * seg) + 1:end) = [];
            id = [id rpls_vec( ones( ng(cs), 1) * (1:seg) )];
            cvpos( cb( cs, id( :, 1)), 1) = id( :, 2);
        end
        %Distribute the remaining samples
        temp = find( cvpos == 0);
        if ~isempty( temp)
            [z, id] = sort( rand( seg, 1) );
            cvpos( temp) = id( 1:length( temp) );%( id) ) = rpls_vec( ones( length( id)/seg, 1) * (1:seg) );
        end
    else
        cvpos=rpls_vec((1:seg)'*ones(1,ceil(r/seg)));
        cvpos=cvpos(1:r);

        cvpos(j)=cvpos; %AAR 270905: This was set wrongly. Corrected it
    end
    
    %AAR 010113 Something was wrong with regards to the replicates
    idt = (1:size( n, 1) )';
    idn = idt( j( :, 1) );
    cvc( idn, 1) = cvpos;
    temp = cvc * ones( 1, size( n, 2) );
    nvec = rpls_vec( n);
    tempvec = rpls_vec( temp);
    tempvec( isnan( nvec) ) = [];
    nvec( isnan( nvec) ) = [];
    ind( nvec) = tempvec;

else
    count=1;
    if length(j)<length(per)
        ind=1;
    else
        for i=1:length(per)
            ind(n(j([count end-(count-1)]),1))=count;
            count=count+1;
        end
        j=j(count:end-(count-1));
        if ~isempty(j)
            seg=seg-(count-1)*2;
            for i=1:length(per)-1
                m=rpls_vec(1:seg);
                id=zeros(seg,1);
                [k,p]=rpls_fjernlike(round(per(i)*m));
                p=p(k(:,1)>0,1);
                ind(n(j(p),1))=i;
                id(p)=1;
                j=j(id==0);
                seg=seg-length(p);
                per=per./sum(per(i+1:end));
            end
            ind(n(j,1))=length(per);
        end
    end
end
% while sum( ind == 0) > 0
%     zero=find(ind==0);
%     ind(zero)=ind(zero-1);
% end