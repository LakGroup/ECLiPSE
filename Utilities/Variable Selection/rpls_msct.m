function [xc,w,xr,xw,b]=rpls_msct(x,ref,ax,pol,wav,w,isc)
%Performs MSC-family spectral corrections
% [xc,w,xr,xw,b]=rpls_msct(x,ref,ax,pol,wav,w,isc)
%
%INPUT:
% x     Original spectra
% ref   Reference spectra
% ax    Axis of spectra
% pol   Degree of reference correction (polynomial)
% wav   Degree of wavelength correction (polynomial)
% w     Weighted version of the algorithm. w = NaN means no weighting
%        (default). A vector can be used for a pre-defined weighting of the
%        wavelengths, while w=-1 means that the algorithm finds the best
%        weights (See ref 4). If set as a positive integer, Loopy MSC is
%        performed
% isc   0 (default) Normal MSC, 1 - Inverse MSC
%
%OUTPUT:
% xc    Corrected spectrum
% w     Weights which were used
% xr    Reference correction
% xw    Wavelength correction
% b     Correction coefficients (Col# 1-3 0-2nd order reference correction,
%        remaining columns are axis corrections (1-nth order)

% Copyright, 2013 - 
% This M-file and the code in it belongs to the holder of the
% copyrights and is made public under the following constraints:
% It must not be changed or modified and code cannot be added.
% The file must be regarded as read-only. 
%
% Åsmund Rinnan
% SPECC
% Department of Food Science
% Faculty of Science
% University of Copenhagen
% Rolighedsvej 30, 1958 Frederiksberg C, Denmark
% e-mail: aar@food.ku.dk

%The MSC part is based on ref 1-4 and the ISC on ref 4-6:
%1. P Geladi, D MacDougal, H Martens: Linearization and scatter correction
%    for near-infrared reflectance spectra of meat, Applied Spectroscopy, 
%    1985, 39, 491-500 
%2. H. Martens and E. Stark: Extenden multiplicative signal correction and
%    spectral interference subtraction: new preprocessing methods for near
%    infrared spectroscopy, Journal of Pharmaceutical and Biomedicinal
%    Analysis, 1991, 9, 625-635
%3. H Martens, JP Nielsen, SB Engelsen: Light scattering and light
%    absorbance separated by extended multiplicative signal correction. The
%    application to near-infrared transmission analysis of powder mixtures,
%    Analytical chemistry, 2003, 75, 394-404 
%4. Neal B. Gallagher, Thomas A. Blake, Paul L. Gassman: Application of
%    extended inverse scatter correction to mid-infrared reflectance
%    spectra of soil, Journal of Chemometrics, 2006, 19 (5-7), 271-281 
%5. DK Pedersen, H Martens, JP Nielsen, SB Engelsen: Near-Infrared
%    Absorption and Scattering Separated by Extended Inverse Signal
%    Correction (EISC): Analysis of Near-Infrared Transmittance Spectra of
%    Single Wheat Seeds, Applied Spectroscopy, 2002, 56, 1206-1214
%6. Helland, Inge S; Næs, Tormod; Isaksson, Tomas: Related versions of the
%    multiplicative scatter correction method for preprocessing 
%    spectroscopic data, Chemometrics and Intelligent Laboratory Systems, 
%    1995 (29): 233-241

%270508 AAR Added the inverse form, although not in the weighted or
%            iterated version
%081007 AAR Some corrections of the quadratic term of the reference
%240907 AAR

[r,c]=size(x);

%Set default values
if nargin<2 || isempty(ref)
    ref=rpls_nanmean(x); %240915 Not optimal
end
if nargin<3 || isempty(ax) || any( isnan( ax) )
    ax=1:size(x,2);
end
if nargin < 4 || isempty( pol)
    pol = 1;
end
if nargin < 5 || isempty( wav)
    wav = 0;
end
if nargin<6 || isempty(w)
    w=NaN;
end
if nargin<7 || isempty(isc)
    isc=0;
end
if length(w)>1 && length(w)~=c
    error('''w'' should either have length 1 or size(x,2)')
end
if length(w)==c
    w=sparse(diag(w));
else %Added to activate LMSC
    if w > 0
        num = w;
    else
        num = 1;
    end    
end

st = 0;
    
while st < num

    if isc==0 % Normal MSC
        X=ones(c,1);
        vid = true( c, 1);
        if pol>0
            X=[X zeros(c,pol)];
            for i=1:pol
                X(:,i+1)=vec(ref).^i;
            end
            vid = vec( ~isnan( ref));
        end
        
        if wav>0
            X=[X zeros(c,wav)];
            for i=1:wav
                X(:,pol+i+1)=vec(ax).^i;
            end
        end
        
        if length(w)==1
            %Finding the corretions
            xi=pinv(X( vid, :),eps*3);
            b=(xi*x( :, vid)')';
            if w == -1
                opt = 1;
            else
                opt = 0;
            end
        elseif length(w)==c
            xi=pinv(X'*w*X,3*eps);
            b=(xi*X'*w*x')';
            opt=0;
        end
        
        %Wavelength correction
        xw=x-b(:,1)*ones(1,c)-b(:,pol+2:end)*X(:,pol+2:end)';
        
        %Reference correction
        switch pol
            case 1
                xc=xw./(b(:,2)*ones(1,c));
            case 2
                %         xc=(-b(:,2)*ones(1,c) + sqrt((b(:,2)*ones(1,c)).^2 - 4*((b(:,1)*ones(1,c)-x).*(b(:,3)*ones(1,c))))) ./ (2*b(:,3)*ones(1,c));
                e4ac = 4 * xw.*( b( :, 3) * ones( 1, c) ); 
                eb = b( :, 2) * ones( 1, c);
                e2a = 2 * b( :, 3) * ones( 1, c);
                xc = (-eb + sqrt( eb.^2 + e4ac) )./ e2a;
                %Just by looking at the spectra it seems that the final correction
                %would be xc=real(xc)-imag(xc); However, I cannot defend this
                %mathematically in any way!
                xc=real(xc)-imag(xc);
            otherwise
                xc=xw;
        end
        
        if opt==1
            %     w=max(var(xc))-var(xc);
            w=1./sum((x-xc).^2);
            w=w/max(w);
            %Should be otpimized :) It actually converges it's possible to
            %check sqrt( mean( (W( :, last) - W( :, previous)).^2) )
            con = 1;
            it = 1;
            while con > 1e-12 && it < 40
                W( :, it) = w;
                w = sparse( diag(w));
                xi = pinv( X'*w*X, 3*eps);
                b = (xi*X'*w*x')';
                xw = x - b( :, 1) * ones( 1, c) - b( :, pol+2:end) * X( :, pol+2:end)';
                if size(b,2)>1
                    xc = xw ./ (b( :, 2) * ones( 1, c));
                else
                    xc = xw;
                end
                w = 1./sum((x-xc).^2);
                %         w=max(var(xc))-var(xc);
                w = w/max(w);
                con = sqrt( mean( ( W( :, end) - w').^2) );
                %             disp( num2str( [ it con] ) )
                it = it + 1;
            end
        end
        
        xr = xc - xw;
        xw = b(:,1) * ones(1,c) + b(:,pol+2:end)*X(:,pol+2:end)';
        
        st = st + 1;
        
    else % Inverse MSC
        X = ones( c, pol+wav+1);
        if wav > 0
            X( :, pol+2:end) = (vec(ax)*ones(1,wav)).^(ones(c,1)*(1:wav));
        end
        for i=1:r
            if pol > 0
                X( :, 2:(pol+1)) = (vec(x(i,:))*ones(1,pol)).^(ones(c,1)*(1:pol));
            end
            xi=pinv(X,eps*3);
            b(i,:)=(xi*ref')'; %Calculating the parameters
            xc(i,:)=X*b(i,:)';  %Correcting the spectra
            id=[1 pol+2:(pol+wav+1)];
            xw(i,:)=X(:,id)*b(i,id)';
            xr(i,:)=xc(i,:)-xw(i,:);
        end
    end
    
    ref = mean( xc);
    st = st + 1;
end