function [ model, evolve] = rpls_NotPLSToolbox( X, Y, opt)
% [ model, evolve] = rpls( X, Y, opt)
%
%INPUT:
% X      X-Matrix, with samples in rows
% Y      Y-matrix, only good for one Y at the moment
% opt    Type opt = rpls to get default settings
%
%OUTPUT:
% model  Includes the model parameters for the model
% evolve Information on how the model has evolved

% Please refer to the paper
% Rinnan, Andersson, Ridder, Engelsen (2014): Recursive weighted partial least squares (rPLS):
% an efficient variable selectionmethod using PLS, Journal of Chemometrics,
% DOI: 10.1002/cem.2582

% Copyright, 2014 - 
% This M-file and the code in it belongs to the holder of the
% copyrights and is made public under the following constraints:
% It must not be changed or modified and code cannot be added.
% The file must be regarded as read-only. 
% In case of doubt, contact the holder of the copyrights.
%
% Åsmund Rinnan
% E-mail asmundrinnan@gmail.com

if nargin == 0
    model.PP.met = {'mc'};
    model.PP.set = NaN;
    model.nLV = 3;
    model.Method = 'rpls';
    model.Val.met = 'CV';
    model.Val.ind = NaN;
    model.Val.rep = NaN;
    model.Val.Info = {'met - Validation method: CV, Boot or Test'; ...
        'ind - Index for CV or Test, or number of Bootstrap samples'; ...
        'rep - Sample number. Replicates have same number'};
    model.Set = [-1 -1];
    model.Info = {'PP - Preprocessing (see ''prepro'')'; 'nLV - Number of factors'; ...
        'Method - rPLS'; 'Ind - Variable index used for iPLS and GA-PLS';...
        'Val - See Val.Info'; 'Set - Two numbers. [Number of iterations, Limit of removing a variable]'; ...
        'Set(1) = Inf: until no more model improvement, default'};
    model.View = false;
    model.Set = [Inf 1e3];
    return
end

opt.Val.met = lower( opt.Val.met);
if strcmp( opt.Val.met, 'cv')
    if isnan( opt.Val.ind)
        opt.Val.ind = 1:size( X, 1); %Make LOO-CV
    end
    cvopt = rpls_cv;
    if opt.nLV == -1
        cvopt.nLV = 30;
    else
        cvopt.nLV = opt.nLV;
    end
    cvopt.xpp.met = opt.PP.met;
    cvopt.xpp.set = opt.PP.set;
    cvopt.view = false;    
end

if opt.Set(2) < 0
    opt.Set(2) = 1e3;
end

[r, c] = size( X);
test = true;
evolve.w = [zeros( 1, c); ones( 1, c)];
evolve.ind = [false( 1, c); true( 1, c)];
evolve.num = [0; c];
evolve.rms = sqrt( mean( rpls_cen_std( Y).^2) );
evolve.nLV = NaN;
% ex = w;
count = 2;
tf = 0;
tr = 1;


while test && count < 100
    Xnew = X.* ( ones( r, 1) * evolve.w( count, :));
    Xcal = X( :, evolve.ind( count, :) );
    Xnew = Xnew( :, evolve.ind( count, :) );
    %Validation method
    switch opt.Val.met
        case 'cv'
            %Make it work with repeated CV
            for cr = 1:size( opt.Val.ind, 2)
                cvopt.Val.ind = opt.Val.ind( :, cr);
                %X is weighted according to B
                mod = rpls_cv( Xnew, Y, cvopt);
                %A variable reduced version of X (no weights)
                modc = rpls_cv( Xcal, Y, cvopt); 
                %Search for the optimal number of LVs
                if opt.nLV == -1
                    rmscal( cr, 1:length( modc.Val.rms( 2:end) ) ) = modc.Val.rms( 2:end);
                    rmstemp( cr, 1:length( mod.Val.rms( 2:end) ) ) = mod.Val.rms( 2:end);
                    Ball{cr} = mod.Val.B;
                else
                    b( :, cr) = rpls_nanmean( mod.Val.B{end}' )'; %Added 091110 AAR
                    rmscal( cr) = modc.Val.rms( end);
                    rmstemp( cr) = mod.Val.rms( end);
                end
            end
            %Search for the optimal number of LVs
            if opt.nLV == -1
                rmstemp( rmstemp == 0) = NaN;
                rmscal( rmscal == 0) = NaN;
                [ ~, evolve.nLV( count)] = min( modc.Val.rms( 2:end) );
                cvopt.nLV = evolve.nLV( count);
%                 evolve.nLV( count + 1) = estfac( modc);
                evolve.detail.rmsall( count + 1, 1:length( rmscal) ) = rmscal;
                evolve.rms( count) = rpls_nanmean( rmscal( :, evolve.nLV( count + 1) ) );
                fac = estfac( mod);
                for cr = 1:size( opt.Val.ind, 2)
                    b( :, cr) = rpls_nanmean( Ball{cr}{fac}, 2);
                end                
                clear rmstemp Ball
            else
                evolve.rms( count) = mean( rmscal);
                evolve.nLV( count) = opt.nLV;
                clear rmstemp
            end
        case 'boot' %Not optimized
            ind = ceil( rand( r, opt.Val.ind) * r);
            for cr = 1:size( ind, 2)
                Xcal = Xnew( ind( :, ci), :);
                Xcal = rpls_prepro( Xcal, opt.PP.met, opt.PP.set);
                Ycal = rpls_cen_std( Y( ind( :, cr) ) );
                
                temp = apls( Xcal, Ycal, opt.nLV);
                b( :, cr) = temp( :, end);
            end
            b = mean( b, 2); %Only extract the average b-values
        case 'test'
        otherwise   
            Xcal = rpls_prepro( Xnew, opt.PP.met, opt.PP.set);
            Ycal = rpls_cen_std( Y);
            
            b = apls( Xcal, Ycal, opt.nLV);
            b = b( :, end);
    end
    
    %Check if there is no improved model performance
    if isinf( opt.Set(1) )
        if evolve.rms( count) > evolve.rms( count - 1)
            tr = tr + 1;
            if tr > 2
                test = false;
            end
        else
            tr = 1;
        end
    end
    
    if test
        count = count + 1;
        b( isnan( b) ) = 0; %Added 091110 AAR
        B = zeros( 1, c);
        B( evolve.ind( count - 1, :)) = b;
        evolve.w( count, :) = evolve.w( count - 1, :).* abs( B);
        evolve.w( count, :) = evolve.w( count, :)/ max( evolve.w( count, :) );
        evolve.ind( count, :) = evolve.w( count, :) > (1/ opt.Set(2));
        evolve.num( count, 1) = sum( evolve.ind( count, :));

        %See if there is no difference in the weights
        temp = sum( abs( diff( evolve.w( count-1:count, :) ) ) );
        if temp < 1e-8 && evolve.num( count - 1, 1) == evolve.num( count, 1) 
            test = false;
        end
        
        %Check if the number of variables do not decrease further
        if test && ( evolve.num( count) == evolve.num( count-1) && evolve.num( count) > 0)
            if c - evolve.num( count) > cvopt.nLV || opt.nLV == -1 || evolve.num( count) == c
                tf = tf + 1;
                %Try to boost the search a bit
                evolve.w( count, :) = evolve.w( count - 1, :).* abs( B).^(tf); 
                evolve.w( count, :) = evolve.w( count, :)/ max( evolve.w( count, :) );
                evolve.ind( count, :) = evolve.w( count, :) > (1/ opt.Set(2));
                evolve.num( count, 1) = sum( evolve.ind( count, :));
%                 if evolve.num( count, 1) < opt.nLV && opt.nLV > 0                    
%                     temp = find( evolve.w( count, :) > 0);
%                     bad = 0;
%                     while length( temp) < opt.nLV
%                         temp = find( evolve.w( count - bad, :) > 0);
%                         bad = bad + 1;
%                     end
%                     [i, j] = sort( evolve.w( count, temp)', 'descend');
%                     evolve.ind( count, :) = 0;
%                     evolve.ind( count, temp( j( 1:opt.nLV) ) ) = 1;
%                     evolve.num( count, 1) = opt.nLV;
% %                     test = false;
%                 end
                if tf > 10
                    evolve.w( end, :) = [];
                    evolve.num( end) = [];
                    evolve.ind( end, :) = [];
                    test = false;
                end
            else
                test = false;
            end
        else
            tf = 0;
        end
        

        if opt.Set(1) > 0 && count > opt.Set(1)
            test = false;
        end
        clear b
    end
    if opt.View
        fprintf( [num2str( evolve.num( end)) ' ' num2str( evolve.rms( end), '%4.3f') char( 10)])
    end
end
[evolve.num, id] = rpls_fjernlike( evolve.num);
evolve.detail.id = id;
evolve.ind = evolve.ind( id( :, 1), :);
evolve.detail.rms = evolve.rms;
evolve.rms = evolve.rms( id( :, 1));
evolve.nLV = evolve.nLV( id( :, 1));
if isinf( opt.Set(1))
    [temp, evolve.optimal] = min( evolve.rms( 2:end), [], 2);
else
    evolve.optimal = size( evolve.ind, 1);
end
Xnew = X( :, evolve.ind( evolve.optimal, :));
for cr = 1:size( opt.Val.ind, 2)
    cvopt.ind = opt.Val.ind( :, cr);
    model(cr) = rpls_cv( Xnew, Y, cvopt);
    model(cr).ind = evolve.ind( evolve.optimal, :);
end