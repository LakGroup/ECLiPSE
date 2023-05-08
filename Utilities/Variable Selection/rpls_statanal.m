function [stat,det] = rpls_statanal(Y,Yp,rep,file,op)
% [stat,det]=rpls_statanal(Y,Yp,rep,file,op);
%
%Calculates the prediction statistics for one variable at a time.
%
% INPUT:
%  Y    Reference values
%  Yp   Predicted values
%  rep  Indicates replicates (see 'nam2rep')
%        Default = (1:length(Y))'
%  file A vector indicating if the Yp are from different files (running
%        numbers)
%        Default = ones(length(Y),1)
%  op   0 if no figures, 1 if figures, 2 if bias adjusted pred vs act
%        Default = 1
%
% OUTPUT:
%  stat Statistics for the prediction
%  det  Optional, gives bias and slope/intercept information in case
%        length(fjernlike(file))>1
%
% See also: apls, cv

%AAR 120912 Added a test to check whether there is a variation in the
%            predicted Y-values (i.e. in order to be able to caclulate the
%            RMSEP-value for the 0th order model)
%AAR 160107 Adjusted the plotting to make it nicer
%AAR 311005 Serious upgrading, dividing the function into several
%            subfunctions
%AAR 271005 Corrected for small files
%AAR 201005 Added file-wise bias and slope/intercept correction
%           Corrected the repeatability
%AAR 060105 Added intercept to the statistics and also bias and intercept
%            adjusted RMSEP values

if min(size(Y))>1
    error('statanal can only handle one reference value at a time')
end

Y=vec(Y);
if size(Yp,1)~=length(Y)
    if size(Yp,2)~=length(Y)
        error('The lengths of Yp and Y are different')
    else
        Yp=Yp';
    end
end

if nargin<3 || isempty(rep)
    rep=(1:length(Y))';    
end
if nargin<4 || isempty(file)
    file=ones(length(Y),1);
end
rep(isnan(Y))=[];
file(isnan(Y))=[];
Yp=Yp(~isnan(Y),:);
Y=Y(~isnan(Y));

if nargin<5 || isempty(op)
    op=1;
end
if length(rep)~=length(Y)
    disp('''rep'' had the wrong length. Set to default')
    rep=(1:length(Y))';
end
if length(file)~=length(Y)
    disp('''file'' had the wrong length. Set to default')
    file=ones(length(Y),1);
end

uf = rpls_fjernlike(file);
Ypred{1}=[];
Ypred{2}=[];
Ypred{3}=[];
Yt=[];
for i=1:length(uf)
    [det(i),Yrep,Yc] = statpred(Y(file==uf(i)),Yp(file==uf(i),:),rep(file==uf(i)));
    for j=1:3
        Ypred{j}=[Ypred{j};Yc{j}];
    end
    Yt=[Yt;Yrep];
    sam(i)=det(i).Samples(1);
end

yvar=size(Yp,2);
stat.RMSEP=sqrt(mean((Yt(:,ones(yvar,1))-Ypred{1}).^2));
%Skewness
% stat.skew=sum((cen_std(Yt(:,ones(yvar,1))).^3)./((std(Ypred{1}).^3)*(size(Ypred{1},1)-1)); 
stat.Bias=mean(Yt(:,ones(yvar,1))-Ypred{1});
%Should check for files with only one sample
stat.SEP=sqrt(sum((Yt(:,ones(yvar,1))-Ypred{2}).^2)/(sum(sam-1)));
for i=1:yvar
    if var( Ypred{1}( :, 1) ) > eps
        temp=polyfit(Yt,Ypred{1}(:,i),1);
        stat.Slope(i)=temp(1);
        stat.Inter(i)=temp(2);
    else
        stat.Slope(i) = NaN;
        stat.Inter(i) = NaN;
    end
end
%Should check for files with less than three samples
stat.SEPCorr=sqrt(sum((Yt(:,ones(yvar,1))-Ypred{3}).^2)/sum(sam-2));
for i=1:3
    temp=corrcoef([Yt Ypred{i}]);
    stat.R2(i,:)=temp(1,2:end).^2;
end
stat.RepAb=repab(Yp,rep);
stat.Samples=[sum(sam) length(Y)];
stat.Range=[min(Y) max(Y)];

stat.R2=full(stat.R2);

if op>0
    for i=1:yvar
        figure        
        plot(Ypred{op}(:,i),Yt,'bo','MarkerSize',4,'MarkerFaceColor','b','LineWidth',2)
        hold on
        plot([floor(min(min(Ypred{op}(:,i)),min(Yt))) ceil(max(max(Ypred{op}(:,i)),max(Yt)))],[floor(min(min(Ypred{op}(:,i)),min(Yt))) ceil(max(max(Ypred{op}(:,i)),max(Yt)))],'k','LineWidth',2)
        xlabel('Predicted')
        ylabel('Actual')
        %Make the axis range a lot nicer
        lim=[min(min(Ypred{op}(:,i)),min(Yt)) max(max(Ypred{op}(:,i),max(Yt)))];
        pow=min(ceil(-log10(lim))+1);
        lim=lim*10^pow;
        lim=[floor(lim(1))/10^pow ceil(lim(2))/10^pow];
        axis([lim lim])
%         axis([floor() ceil() floor() ceil(max(max(Ypred{op}(:,i)),max(Yt)))])

        %Plot the residual AFTER correction
        res=Yt-Ypred{op}(:,i);

        stp=rangeX(res)/25;
        bin=sort([(stp/2:stp:ceil(max(res)))';(-stp/2:-stp:floor(min(res)))']);
        bin=[bin zeros(size(bin,1),1)];
        for j=1:length(bin)-1
            bin(j,2)=length(find(res>=bin(j,1) & res<bin(j+1,1)));
        end
        a=1;
        while bin(a,2)==0
            a=a+1;
        end
        b=size(bin,1);
        while bin(b,2)==0
            b=b-1;
        end
        figure
        h=bar(bin(a:b,1)+stp/2,bin(a:b,2));
        set(h,'FaceColor',[1 1 1]) % Has to be phased out for JJ
        % set(gca,'XGrid','on') % Also wanted by JJ
        % set(gca,'YGrid','on') % Also wanted by JJ
        xlabel('Actual - Predicted')
        ylabel('Number of samples')
        x=bin(a:b,1)+stp/2;
        x=min(x):rangeX(x)/1000:max(x);
        %In case there are replicates this has to be accounted for
        d=npdf(x,mean(res),std(res));

        hold on;plot(x,(d./max(d))*max(bin(:,2)),'r','LineWidth',2)
        %         hold on;plot(x,(d./max(d))*max(bin(:,2)),'r','LineWidth',1) %JJ wants linewidht = 1
    end
end

%---------------------------------------
function xnew=vec(x)
%xnew=vec(x)
%
% Vectorices the matrix x

xnew=x(:);

%---------------------------------------
function [part,Yrep,Yc]=statpred(Y,Yp,rep)
%[part,Yrep,Yc]=statpred(Y,Yp,rep)
%
% INPUT
%  Y     Measured values
%  Yp    Predicted values
%  rep   Vector with same number for replicates
%
% OUTPUT
%  part  Struct with RMSEP, SEP, SEPCorr, Slope, Intercept, Bias, R^2,
%         Samples, Range and RepAb
%  Yrep  Reference Y. Replicates removed (if present)
%  Yc    Cell-array. (1) Uncorrected
%                    (2) Bias corrected
%                    (3) Slope corrected

[Ym,Yr]=rpls_repmean([Y Yp],rep);
p=size(Ym,2)-1;
Yrep=Ym(:,1);
part.RMSEP=sqrt(mean((Ym(:,ones(p,1))-Ym(:,2:end)).^2));
% part.skew=sum(cen_std(Ym(:,2:end)).^3)./(std(Ym(:,2:end)).^3*(size(Ym,1)-1));
% %Skewness
Yc{1}=Ym(:,2:end);
temp=corrcoef(Ym);
part.R2=temp(1,2:end);
if size(Ym,1)>2
    part.Bias=mean(Ym(:,ones(p,1))-Ym(:,2:end));
    Yc{2}=Ym(:,2:end)+part.Bias(ones(size(Ym,1),1),:);
    part.SEP=std(Ym(:,ones(p,1))-Ym(:,2:end));
else
    part.Bias=NaN;
    part.SEP=NaN;
    Yc{2}=Ym(:,2:end);
end
if size(Ym,1)>3
    for i=1:p
        if var( Ym( :, i + 1) ) > eps
            %181012 AAR Probably need to do a check-up here on the dimensionality
            temp=polyfit(Ym(:,i+1),Ym(:,1),1);
            part.Slope(i)=temp(1);
            part.Inter(i)=temp(2);
            Ym(:,i+1)=Ym(:,i+1)*temp(1)+temp(2);
        else
            part.Slope(i) = NaN;
            part.Inter(i) = NaN;
            Ym( :, i+1) = Ym( :, i+1) * NaN;
        end
    end
    part.SEPCorr=sqrt(sum((Ym(:,ones(p,1))-Ym(:,2:end)).^2)/(size(Ym,1)-2));
else
    part.Slope(1:p)=NaN;
    part.Inter(1:p)=NaN;
    part.SEPCorr=NaN;
end
Yc{3}=Ym(:,2:end);
part.RepAb=repab(Yp,rep);
part.Samples=rpls_fjernlike([size(Ym,1);size(Y,1)])';
part.Range=[min(Y) max(Y)];

%---------------------------------------------
function x=repab(Yp,rep)

[temp,Yr]=rpls_repmean(Yp,rep);
if size(Yp,1)~=size(temp,1)
    [j,k]=rpls_fjernlike(rep);
    %Remove those samples that only have one replicate
    ind=find(~isnan(k(:,2)));
    ind=vec(k(ind,:)');
    ind=ind(~isnan(ind));
    %Calculate the difference from mean
    res=(Yr(ind,:)-Yp(ind,:)).^2;
    %Find replicates
    [f,g]=rpls_fjernlike(rep(ind));
    %In case there is a mixture of samples with and without replicates,
    %the samples without replicates have to be removed
    f(isnan(g(:,2)))=[];
    g(isnan(g(:,2)),:)=[];
    %Set index for sparse matrix
    u=vec((1:length(f))'*ones(1,size(g,2)));
    t=vec(ones(size(f))*(1:size(g,2)));
    f=vec(g);
    u(isnan(f))=[];
    t(isnan(f))=[];
    f(isnan(f))=[];
    %Calculate repeatability per column in Yp
    for i=1:size(Yp,2)
        w=sparse(u,t,res(f,i));
        w=sqrt(sum(w')'./(sum(~isnan(g'))'-1));
        x(i)=full(mean(w));
    end
else
    x=ones(1,size(Yp,2))*NaN;
end