function model = rpls_cv( X, Y, opt)
%model = cv( X, Y, opt)
% 
% Cross validation by pls, pcr or PARAFAC. Leave-one-out or
% other cross validation models may be used.
%
% INPUT
%  X        X-array
%  Y        Response variables
%  opt      Includes the fields: nLV, met, ind, pp, fac and view. 
%            Type 'opt = cv;' to get the default values
%
% OUTPUT
%  model    A struct array with two main fields: 'Val' and 'Cal'. These
%            again contains the following fields: calibration model,
%            predictions, residuals in X
%
% See alos: cvind

% Copyright, 2014 - 
% This M-file and the code in it belongs to the holder of the
% copyrights and is made public under the following constraints:
% It must not be changed or modified and code cannot be added.
% The file must be regarded as read-only. 
% In case of doubt, contact the holder of the copyrights.
%
% Åsmund Rinnan
% E-mail asmundrinnan@gmail.com

if nargin == 0 || nargin < 3
    model.nLV = 10;
    model.class = false;
    model.method = 'pls';
    model.weights.X = NaN;
    model.weights.Y = NaN;
    model.Val.ind = NaN;
    model.Val.rep = NaN;
    model.xpp.met{1} = 'mc';
    model.xpp.set(1).ref = NaN;
    model.xpp.set(1).ax = NaN;
    model.xpp.set(1).set = NaN;
    model.ypp.met{1} = 'mc';
    model.ypp.set(1).ref = NaN;
    model.ypp.set(1).ax = NaN;
    model.ypp.set(1).set = NaN;    
    model.fac = NaN;
    model.view = 1;
    model.info = {'nLV - Number of Latent Variables to be used'; ...
        'class - false for continous data, true for classification problems'; ...
        'method can be ''PLS'', ''WPLS'', ''PARAFAC'', ''PCR'' or ''N-PLS'''; ...
        'Val - ind: Index for cross-validation (see ''cvind'' or cvpos''';...
        '      rep: vector of sample numbers (i.e. replicates have same number';...
        'xpp - Pre-processing of X. It will always be done in the following order:'; ...
        '       SNV/ MSC, followed by SG, followed by scaling and centering'; ...
        '       If multiple pp is wanted, give them in cells.'; ...
        'ypp - As above but for Y.'; ...
        'fac - Initial guesses of the PARAFAC factors.'; ...
        'view - 1 = show progress, anything else = don''t show progress'};
    if nargin == 0
        return
    end
end

%Check if X and Y has same number of rows
n = size(X);
[rY,kY]=size(Y);
if n(1)~=rY
    error('''X'' and ''Y'' should have the same number of samples.')
end

if any( isnan( Y) )
    error( 'You have missing values in y, remove these prior to running ''cv''')
end

%Reduce the number of max factors in case size(X) is small
if opt.nLV>min(size(X))
    opt.nLV = min(size(X));
end

%Check if the method is valid
opt.method = lower( opt.method);
if sum( strcmp( opt.method, {'pcr', 'pls', 'parafac', 'npls', 'wpls'}) ) == 0
    error( [opt.method ' is not a recognized method. Select between: pcr, pls, parafac, npls and wpls'])
end

if isnan( opt.Val.rep)
    opt.Val.rep = ( 1:n(1) )';
end

%The default value of 'pos'
if isnan( opt.Val.ind)
    opt.Val.ind = cvind( y, opt.Val.rep, min( n(1),20), true( 1, 2));
end

%In case of PLS-DA and the Y-variable is just given as a classindex rather
%than a dummy variable
if opt.class && size( Y, 2) == 1 && length( rpls_fjernlike( Y( :, 1) ) ) > 2    
    [~, yi] = rpls_fjernlike( Y);
    yc = (1:size( yi, 1) )' * ones( 1, size( yi, 2) );
    yc( isnan( yi) ) = NaN;
    ytemp = zeros( size( Y) );
    yi = vec( yi');
    yc = vec( yc');
    yi( isnan( yi) ) = [];
    yc( isnan( yc) ) = [];
    ytemp( yi) = yc;
    Y = full( sparse( (1:size( Y, 1) )', ytemp, true( size( Y, 1), 1) ) );
end

%If fac is not given, and method is PARAFAC, this is calculated
if strcmp( 'par', opt.method)
   if isnan( opt.fac)
      opt.fac = parafac(X,nLV,[0 0 0 0 -1]);
   end
end

opt.view = opt.view == 1;

%There may be columns in X or Y with zero variance
if any( var( Y) == 0) || any( var( X) == 0)
    error( 'Some of your variables have zero variance! Remove them, and run again')
end

%Initialize the prediction of Y
ypred = zeros( rY, opt.nLV);

%Start the waiting bar
if opt.view
    h = waitbar(0,['Cross validation is now running using ' opt.method]);
    set(h,'units','normalized');
    set(h,'position',[.4 .6 .3 .075]);
end

%Initialize the output struct
model.Cal=struct('Yp',[],'rms',[],'Xres',[],'B',[],'T',[],'P',[],'U',[],'C',[],'W',[]);
model.Val=struct('Yp',[],'rms',[],'Xres',[],'B',[],'T',[], 'U', []);
model.Org.Y = Y;
model.Org.rep = opt.Val.rep;

%Pre-processing order
for pp = {'xpp', 'ypp'}
    temp = rpls_vec( opt.(pp{1}).met);
    ind = zeros( 1, 3);    
    ind( 1) = max( [0; find( strcmp( temp, 'snv') + strcmp( temp, 'msc') )] );
    ind( 2) = max( [0; find( strcmp( temp, 'sg') )] );
    ind( 3) = max( [0; find( strcmp( temp, 'mc') + strcmp( temp, 'as') + strcmp( temp, 'ps') )] );
    ind( ind == 0) = [];
    opt.(pp{1}).met = opt.(pp{1}).met( ind);
    opt.(pp{1}).set = opt.(pp{1}).set( ind);
end

%Indicator if there are missing values in the original X-matrix
test = isnan( sum( rpls_vec( X) ) ); 

%The cross-validation itself
ur = rpls_fjernlike( opt.Val.ind);
for ci = 1:length( ur)
    if opt.view == 1
        waitbar( ci/length(ur),h)
    end
    Xc = reshape( X( opt.Val.ind ~= ur(ci), :), [sum( opt.Val.ind ~= ur(ci) ) n(2:end)]);
    Xp = reshape( X( opt.Val.ind == ur(ci), :), [sum( opt.Val.ind == ur( ci) ) n(2:end)]);
    Yc = Y( opt.Val.ind ~= ur( ci), :);
    Yp = Y( opt.Val.ind == ur( ci), :);
        
    %The part model and prediction is performed in a separate function
    mod = rpls_regres( Xc, Xp, Yc, Yp, opt);
        
    switch opt.method
        case 'pls'
            if isempty( mod) %Check if any model has been made
                error( 'There are several samples where the variable(s) do not vary')
            else
                %Check if this is a classification problem, and transform the predictions
                %into classifications
                if opt.class                    
                    if size( Y, 2) > 1 %Multi-class problem
                        ypclass = zeros( size( mod.yp{1}) );
                        for cf = 1:size( mod.yp{1}, 2)
                            ypred = zeros( size( mod.yp{1}, 1), length( mod.yp) );
                            for cc = 1:length( mod.yp)
                                ypred( :, cc) = mod.yp{cc}( :, cf);
                            end
                            [~, j] = min( abs( 1 - ypred), [], 2 ); %Find the column closest to one
%030212 AAR Doesn't seem like it makes sense to set some to non-classified.
%            Tested on fish data
%                             ypred = sort( ypred, 2);
%                             j( ypred( :, end) - ypred( :, end-1) < .025) = NaN;
                            ypclass( :, cf) = j;
                        end
                    else %Two class problem
%                         error( 'Should calculate a ROC-curve and thus get a better result!')
                        ypclass = mod.yp > sum( Yc == 0)/ length( Yc);
                        %300812 AAR Changed the classification limit from
                        %            0.5 to be dependent on the number of
                        %            samples in each group
                    end
                    mod.yp = ypclass;
                end
                
                if ~isnan( sum( rpls_vec( Xc) ) ) && test && ~iscell( mod.Tp)
                    mod.Tp = mat2cell( mod.Tp, size( mod.Tp, 1), ones( size( mod.Tp, 2), 1) );
                    mod.Up = mat2cell( mod.Up, size( mod.Up, 1), ones( size( mod.Up, 2), 1) );
                end
                if test
                    model.Val.Yp( opt.Val.ind == ur( ci), 1:size( mod.yp, 2) ) = mod.yp;
                    for cf = 1:size( mod.B, 2)
                        %AAR 091110 Two lines added
                        model.Val.B{cf}( mod.remid == 0, ci) = mod.B( :, cf);
                        model.Val.B{cf}( mod.remid == 1, ci) = NaN;
                    end
                    model.Val.Xres( opt.Val.ind == ur( ci), 1:size( mod.Xres, 2) ) = mod.Xres;
                    if iscell( mod.Tp)
                        for cf = 1:min( opt.nLV, size( mod.Tp, 2) )
                            %Sometimes the number of extracted factors don't match
                            if cf == size( mod.Tp{ cf}, 2)
                                model.Val.T{cf}( opt.Val.ind == ur(ci), 1:cf) = mod.Tp{ cf};
                                model.Val.U{cf}( opt.Val.ind == ur(ci), 1:cf) = mod.Up{ cf};
                            else
                                opt.nLV = cf - 1;
                            end
                        end
                    else
                        for cf = 1:min( opt.nLV, size( mod.Tp, 2) )
                            model.Val.T{cf}( opt.Val.ind == ur( ci), :) = mod.Tp( :, 1:cf);
                            model.Val.U{cf}( opt.Val.ind == ur( ci), :) = mod.Up( :, 1:cf);
                        end
                    end
                else
                    if iscell( mod.yp)
                        Yid = find( ~mod.Yrem);
                        for cy = 1:size( mod.yp, 2)
                            model.Val.Yp{ Yid( cy)}( opt.Val.ind == ur( ci), 1:size( mod.yp{cy}, 2) ) = mod.yp{cy};
                        end
                    else
                        model.Val.Yp( opt.Val.ind == ur( ci), 1:size(mod.yp,2) ) = mod.yp;
                    end
                    if iscell( mod.B)
                        for cy = 1:length( mod.B)
                            for cf = 1:size( mod.B{ cy}, 2)
                                model.Val.B{ cy, cf}( ~mod.remid, ci) = mod.B{ cy}( :, cf);
                                model.Val.B{ cy, cf}( mod.remid, ci) = NaN;
                            end
                        end
                    else
                        for cf = 1:size( mod.B, 2);
                            model.Val.B{cf}( ~mod.remid, ci) = mod.B( :, cf);
                            model.Val.B{cf}( mod.remid, ci) = NaN;
                        end
                    end
                    %Removed the following if statement, and the for loop
                    %inside the 'if'
%                     if size( mod.Tp, 2) ~= opt.nLV%) > numel( mod.B) %I don't know what this test is about :(
%                         for cf = 1:size( mod.Tp, 2)
%                             model.Val.T{cf}( opt.Val.ind == ur( ci), 1:cf) = mod.Tp( :, 1:cf);
%                             model.Val.U{cf}( opt.Val.ind == ur( ci), 1:cf) = mod.Up( :, 1:cf);
%                         end
%                     else
                        model.Val.T( opt.Val.ind == ur(ci), 1:size(mod.Tp, 2) )= mod.Tp; %271110 AAR Changed to cell reference
                        model.Val.U( opt.Val.ind == ur(ci), 1:size( mod.Up, 2) ) = mod.Up; %271110 AAR Changed to cell reference
%                     end
                    model.Val.Xres( opt.Val.ind == ur(ci), 1:size(mod.Xres,2) ) = mod.Xres;
                end
            end
        case 'npls'
    end
            
end

%Close the waitbar
if opt.view
    close(h)
end

%Calculate the prediction error of the validation
if iscell( model.Val.Yp)
    for cy = 1:length( model.Val.Yp)
        model.Val.Yp{ cy}( model.Val.Yp{ cy} == 0) = NaN;
%         id( cy, :) = isnan( sum( model.Val.Yp{ cy}) );
%         model.Val.Yp{ cy}( :, id( cy, :) ) = [];
        factors( cy) = size( model.Val.Yp{ cy}, 2);
        if factors( cy) > 0
            stat = statanal( model.Org.Y( :, cy), model.Val.Yp{ cy}, model.Org.rep, [], 0);
            model.Val.rms( cy, 2:factors( cy) + 1) = stat.RMSEP;
        else
            model.Val.rms( cy, :) = NaN;
        end
    end
%     id = sum( id) > 1;
else
%     if opt.class %Classification problem
%         model.Val.Yp( model.Val.Yp == 0) = NaN;
%         id = isnan( sum( model.Val.Yp) );
%         model.Val.Yp( :, id) = [];
%     else
%         id = false( size( model.Val.Yp, 2), 1);
%     end
    %Find the number of factors for the model
    factors = size( model.Val.Yp, 2);
    
    if opt.class
        if size( Y, 2) > 1 %Multi-class problem
            %Converting the predictions into dummy variables and
            %calculating the misclassification error
            for cf = 1:factors
                dumy = false( size( model.Val.Yp, 1), size( Y, 2) );
                for cc = 1:size( Y, 2)
                    dumy( model.Val.Yp( :, cf) == cc, cc) = true;
                end
                %Should take into account the replicates here as well...
                model.Val.rms( cf + 1) = sum( sum( abs( dumy - Y), 2) ~= 0);
            end
        else %Two-class problem
            model.Val.rms( 2:factors + 1) = sum( model.Val.Yp - Y( :, ones( 1, factors) ) ~= 0);
        end
    else
        stat = rpls_statanal( model.Org.Y, model.Val.Yp, model.Org.rep, [], 0);
        model.Val.rms(2:factors+1) = stat.RMSEP;
    end
end
% if iscell( model.Val.T)
%     model.Val.T = model.Val.T( id( 1:length( model.Val.T) ) == 0 );
% % else
% %     model.Val.T( :, id) = [];
% end
% model.Val.Xres( :, id) = [];

%Making the model for the complete data - only works for PLS
%--------------------------------------
if strcmp( opt.method, 'pls')   
    %Preprocessing
    %X-matrix
    Xt = X;
    for cp = 1:length( opt.xpp.met)
        Xt = rpls_prepro( Xt, opt.xpp.met{cp}, opt.xpp.set(cp));
    end
    %Y-matrix
    Yt = Y;
    for cp = 1:length( opt.ypp.met)
        [Yt, Yps] = rpls_prepro( Yt, opt.ypp.met, opt.ypp.set);
    end
    
    if opt.class
        %110413 AAR New suggestion
        model.Cal.rms = length( rpls_fjernlike( model.Org.rep) )/ size( Y, 2) * (size( Y, 2) - 1);
%         model.Cal.rms = length( rpls_fjernlike( model.Org.rep) );
        model.Val.rms(1) = model.Cal.rms;
        clear ypclass
    else
        %RMSE should be calculated based on the original and NOT
        %pre-processed Y-values
        model.Cal.rms(:, 1) = sqrt( mean( rpls_cen_std( rpls_repmean( Y, model.Org.rep) ).^2) );
        model.Val.rms(:, 1) = model.Cal.rms;
    end
    
    if iscell( model.Val.Yp) %Not done yet, decided to run several PLS1 instead afterall...
        [model.Cal.B, model.Cal.W, model.Cal.T, model.Cal.P, model.Cal.U, model.Cal.C] ...
            = apls( Xt, Yt, opt.nLV, [], [], 1);
        [model.Cal.Yp, ~, model.Cal.Xres] = rpls_plspred( Xt, model.Cal.B, model.Cal.P, Yps.ref);
        for cy = 1:length( model.Cal.Yp)
            stat = rpls_statanal( model.Org.Y( :, cy), model.Cal.Yp{cy}, model.Org.rep, [], 0);
            model.Cal.rms( cy, 2:length( stat.RMSEP) + 1) = stat.RMSEP;
        end
    else
        nLV = size( model.Val.Yp, 2); % In case max nLV was reached for one or more segments
        if length( Yps.ref) == 1
            Yps.ref( 2) = 1;
        end
        if sum(rpls_vec(isnan(Xt)))>0
            for cf = 1:opt.nLV
                mod = rpls_plsmiss( Xt, Yt, cf);
                model.Cal.B( :, cf) = mod.B( end, :)';
                model.Cal.W{cf} = mod.W;
                model.Cal.T{cf} = mod.Tx;
                model.Cal.P{cf} = mod.Px;
                model.Cal.U{cf} = mod.Ty;
                model.Cal.C{cf} = mod.Py;
%                 [model.Cal.B( :, cf), W, T, P, U, C] = plsmiss( Xt, Yt, cf);
%                 model.Cal.W{cf} = W;
%                 model.Cal.T{cf} = T;
%                 model.Cal.P{cf} = P;
%                 model.Cal.U{cf} = U;
%                 model.Cal.C{cf} = C;
                %To find the rms0-value to input into plspred
                temp = Xt;
                tp = mod.Tx * mod.Px';
                temp( isnan( temp) ) = tp( isnan( temp) );
                model.Cal.Yp( :, cf) = temp * model.Cal.B( :, cf) * Yps.ref(2) + Yps.ref(1);
                temp = Xt - mod.Tx * mod.Px';
                temp( isnan( Xt) ) = 0;                
                model.Cal.Xres( :, cf) = sum( temp'.^2);
%                 [model.Cal.Yp( :, cf), T, model.Cal.Xres( :, cf)]...
%                     = plspred( Xt, model.Cal.B( :, cf), model.Cal.P{cf}, Yps.ref);
            end
        else
            if size( Yt, 2) == 1
                [model.Cal.B,model.Cal.W,model.Cal.T,model.Cal.P,model.Cal.U,model.Cal.C] ...
                    = rpls_bipls( Xt, Yt, opt.nLV);
            else
                mod = rpls_simpls( Xt, Yt, opt.nLV);
                model.Cal.B = mod.B;
                model.Cal.T = mod.Tx;
                model.Cal.P = mod.Px;
                model.Cal.W = mod.W;
                model.Cal.U = mod.Ty;
                model.Cal.C = mod.Py;
            end
            [model.Cal.Yp, T, model.Cal.Xres] = rpls_plspred( Xt, model.Cal.B, model.Cal.P, Yps.ref);
        end
        factors = size(model.Cal.T,2);
        if opt.class
            if size( Y, 2) > 1 %Multi-class problem
                for cf = 1:size( model.Cal.Yp{1}, 2)
                    ypred = zeros( size( model.Cal.Yp{1}, 1), length( model.Cal.Yp) );
                    for cc = 1:length( model.Cal.Yp)
                        ypred( :, cc) = model.Cal.Yp{cc}( :, cf);
                    end
                    [i, j] = min( abs( 1 - ypred), [], 2 ); %Find the column closest to one
                    ypclass( :, cf) = j;
                end
            else %Two class problem
                ypclass = model.Cal.Yp > .5;
            end
            model.Cal.Yp = ypclass;
            if size( Y, 2) > 1 %Multi-class problem
                %Converting the predictions into dummy variables and
                %calculating the misclassification error
                for cf = 1:factors
                    dumy = false( size( model.Cal.Yp, 1), size( Y, 2) );
                    for cc = 1:size( Y, 2)
                        dumy( model.Cal.Yp( :, cf) == cc, cc) = true;
                    end
                    model.Cal.rms( cf + 1) = sum( sum( abs( dumy - Y), 2) ~= 0);
                end
            else %Two-class problem
                model.Cal.rms( 2:factors + 1) = sum( model.Cal.Yp - Y( :, ones( 1, factors) ) ~= 0);
            end
            model.Cal.rms(1) = size( Y, 1) * .5;
            model.Val.rms(1) = model.Cal.rms(1);
        else
            stat = rpls_statanal( model.Org.Y, model.Cal.Yp, model.Org.rep, [], 0);
            model.Cal.rms( 2:factors+1) = stat.RMSEP;
        end
    end
end

%I would like to add the Hotelling T2 values
for cf = 1:factors
    if iscell( model.Cal.T)
        model.Cal.H( :, cf) = diag( model.Cal.T{cf} * inv( model.Cal.T{cf}' * model.Cal.T{cf}) * model.Cal.T{cf}');
        try
            model.Val.H( : ,cf) = diag( model.Val.T{cf} * inv( model.Cal.T{cf}' * model.Cal.T{cf}) * model.Val.T{cf}');
        catch
            model.Val.H( :, cf) = NaN;
        end
    else
        model.Cal.H( :, cf) = diag( model.Cal.T( :, 1:cf) * inv( model.Cal.T( :, 1:cf)' * model.Cal.T( :, 1:cf) ) * model.Cal.T( :, 1:cf)' );
        try
            model.Val.H( :, cf) = diag( model.Val.T( :, 1:cf) * inv( model.Val.T( :, 1:cf)' * model.Val.T( :, 1:cf) ) * model.Val.T( :, 1:cf)' );
        catch
            model.Val.H( :, cf) = NaN;
        end
    end
end
% temp = expvar( Xt, model.Cal.T, model.Cal.P);
% model.EV( :, 1) = temp( :, 2);
% temp = expvar( Xt, model.Val.T, model.Cal.P);
% model.EV( :, 2) = temp( :, 2);

model.Org.opt = opt;