function [Xnew,Xm,Xs] = rpls_cen_std(X,m,s)
%function [Xnew,Xm,Xs] = rpls_cen_std(X,m,s)
% 
% This function both mean centers and standardizes a dataset, either
% based on itself, or based on some previously given values.
%
% INPUT:
%  X     The original matrix
%  m     If this is given, these mean values are used
%         Can optionally be set to 'a' (autoscaling), 'c' (mean-centering, Default)
%         or 'p' (Paareto scaling)
%  s     If this is given, the matrix is standardized using
%         these values
%
% OUTPUT:
%  Xnew  The centered/ standardized matrix
%  Xm    The mean values are given here
%  Xs    If asked for the matrix is standardized, and the
%         standard deviations are given here

% 140312 AAR Can now also do Paareto-scaling
% 160508 AAR Can now run without all outputs if autoscaling is wanted
% 080307 AAR Minor changes to take into account if there are none or only 1
%             value for a specific variable
% 010307 AAR Changed the scaling a bit so that it removes those variables
%             with 0 in standard deviation

%Try to rewrite the whole function
[r, c] = size( X);

Xs = ones( 1, c);
opt = 0;
if exist('m', 'var')
    if ischar(m)
        opt=find(strcmp({'c';'a';'p'},m));
        if isempty(opt) % Mean-centering is the default
            opt = 1;
        end
        Xm = rpls_nanmean( X);
        if opt > 1
            Xs = rpls_nanstd( X);
            if opt == 3
                Xs = sqrt( Xs);
            end
        end
    else
        Xm = m;
        if nargin == 3
            Xs = s;
        end
    end
else
    Xm = rpls_nanmean( X);
end

Xnew = X - Xm( ones( r, 1), :);
Xnew = Xnew./ Xs( ones( r, 1), :);
% Xm = rpls_nanmean( X);
% 
% Xs = rpls_nanstd( X);
% 
% 
% [r,c]=size(X);
% E=~isnan(X);
% 
% %In case there are only missing values, these are not changed
% sE=sum(E);
% sE(sE==0)=NaN;
% 
% 
% if opt>0 || nargin==1 || nargout>1
%     Xr=zeros(r,c);
%     Xr(E)=X(E);
%     Xm=sum(Xr)./sE;
%     Xnew=X-ones(r,1)*Xm;
%     if opt > 1 || nargout == 3
%         Xr = zeros( r, c);
%         Xr( E) = Xnew( E);
%         sE(sE==1)=NaN; %If there's only one value, autoscaling will remove the variable
%         Xs=sqrt(sum(Xr.^2)./(sE-1));        
%         Xs(Xs==0)=NaN; %Remove the variables with 0 in std
%         if opt == 3
%             Xs = sqrt( Xs);
%         end
%         Xnew=Xnew(:,~isnan(Xs))./(ones(r,1)*Xs(~isnan(Xs)));
%     end
% else
%     Xnew=X-ones(r,1)*m;
%     if nargin==3        
%         Xnew=Xnew(:,~isnan(s))./(ones(r,1)*s(~isnan(s)));
%     end
% end