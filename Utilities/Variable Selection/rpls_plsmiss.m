function model = rpls_plsmiss( X, Y, nLV, opt)
% model = rpls_plsmiss( X, Y, nLV, opt)
%PLS with missing data
%
%INPUT:
% X      The X-data, they should NOT be preprocessed, see 'opt'
% Y      The Y-data, should not be pre-processed, see 'opt'
% nLV    The number of latent variables to use. Remember that this affects
%         the imputations of the missing values in X
% opt    Type 'opt = rpls_plsmiss' to get default values
%
%OUTPUT:
% model  The resulting PLS model, and the imputations of X
%
% See also: pcamiss, cv, simpls, bipls

%090320 AAR Seriously updated this function including correct
%            pre-processing
%240713 AAR Updated this function so that it runs using 'simpls' rather than 
%            'apls' 

%The default values
if nargin == 0
    model.xpp.met = { 'mc'};
    model.xpp.set = prepro;
    model.xpp.id = NaN;
    model.ypp.met = { 'mc'};
    model.ypp.set = prepro;
    model.tol = [ 1e-9 1e5];
    model.view = false;
    model.Info = { 'xpp - See ''prepro'' for further information';
        'ypp - See ''prepro'' for further information, only needed if you haven''t pre-processed';
        'tol - Two number: The tolerance limit for the change in the iterations, and the max # of iterations';
        'view - See (true) or not (false) the progress'};
    return
end

%In case the default value is missing
if ~isfield( opt.xpp, 'id')
    opt.xpp.id = NaN;
end

if nargin < 3
    opt = rpls_plsmiss;
end


%Finding the missing values in X
Cx = isnan(X);
Cy = isnan( Y);
if any( sum( Cx, 2) == size( Cx, 2) | sum( Cy, 2) == size( Cy, 2) )
    disp( 'Some samples have been removed due to only missing values')
    id = sum( Cx, 2) == size( Cx, 2) | sum( Cy, 2) == size( Cy, 2);
    X( id, :) = [];
    Y( id, :) = [];
    C( id, :) = [];
    model.id = id;
end

[rx, cx] = size(X);
[ ry, cy] = size( Y);
Cx = isnan(X);
Dx = ~Cx;
Cy = isnan( Y);
Dy = ~Cy;


% First guess on the missing values
%X
cxv = ones( rx, 1) * rpls_nanvar( X);
rxv = rpls_nanvar(X')' * ones(1, cx);
temp = cxv - rxv;
mv = ones( rx, 1) * rpls_nanmean(X);
mr = rpls_nanmean(X')' * ones(1, cx);
iX = find( Cx );
id = temp( iX ) > 0;
if isnan( opt.xpp.id)
    X( iX(id) ) = mr( iX(id) );
    X( iX(~id) ) = mv( iX(~id) );
else
    X( iX) = 0;
end
xp = ones( size(X) );
%Y
if sum( isnan( Y)) > 0
    cyv = ones( ry, 1) * rpls_nanvar( Y);
    ryv = rpls_nanvar( Y')' * ones( 1, cy);
    temp = cyv - ryv;
    mv = ones( ry, 1) * rpls_nanmean( Y);
    mr = rpls_nanmean( Y')' * ones( 1, cy);
    iY = find( Cy);
    idy = temp( iY) > 0;
    Y( iY( idy) ) = mr( iY( idy) );
    Y( iX( ~idy) ) = mv( iY( ~idy) );
    missY = true;
    yp = ones( size( Y) );
else
    missY = false;
end

it = 1;
if missY
    diff = mean( (xp(iX) - X(iX)).^2 ) + mean( ( yp( iY) - Y( iY) ).^2);
else
    diff = max( mean( (xp(iX) - X(iX)).^2 ), opt.tol(1) * 2);
end

test = true;
while diff > opt.tol(1) && it < opt.tol(2) && test
    xc = X;
    yc = Y;

    for cm = 1:length( opt.xpp.met)
        [ xc, xppset{cm}] = rpls_prepro( xc, opt.xpp.met{cm}, opt.xpp.set(cm) );
    end
    for cm = 1:length( opt.ypp.met)
        [ yc, yppset{cm}] = rpls_prepro( yc, opt.ypp.met{cm}, opt.ypp.set(cm) );
    end
    
    if size( Y, 2) == 1
        model = rpls_bipls( xc, yc, nLV);
    else
        model = rpls_simpls( xc, yc, nLV);
    end
    
    try
        xp = revpp( model.Tx, model.Px, opt.xpp.met, xppset);
    catch
        xp = model.Tx * model.Px';
    end
    
    try
        yp = revpp( model.Tx, model.Py, opt.ypp.met, yppset);
    catch
        yp = model.Tx * model.Py';
    end
    
    res = X - xp;

    if missY
        diff = mean( (xp(iX) - X(iX)).^2 ) + mean( ( yp( iY) - Y( iY) ).^2);
        Y( iY) = yp( iY);
    else
        diff = mean( (xp(iX) - X(iX)).^2 );
        d( it) = diff;
    end
    X( iX ) = xp( iX );
    
    if opt.view && ceil( it/1e4) == it/1e4
        fprintf( [num2str(it) ': ' num2str( diff, '%3.3e') '\n'] )
    end
    it = it + 1;
end

if diff >= opt.tol(1) || it >= opt.tol(2)
    disp('The solution is only suboptimal, no convergence reached, increase the tolerance limits')
end

model.Xnew = X;
model.iter = it;
model.xres = sum( vec( xc).^2);
%Sometimes nLV is too large ...
nLV = size( model.Tx, 2);
for cf = 1:nLV
    model.xres( cf + 1) = sum( vec( xc - model.Tx( :, 1:cf) * model.Px( :, 1:cf)').^2);
end
model.Xres = sum( ( xc - model.Tx * model.Px').^2, 2);
model.Yp = yp;

try
    model.xppset = xppset;
    model.yppset = yppset;
catch
    model.xppset = NaN;
    model.yppset = NaN;
end