function model = rpls_regres( Xcal, Xval, Ycal, Yval, opt)
% Performs a regression analysis based on a calibration set.
%
%model=regres(Xcal, Xval, Ycal, Yval, opt);
%
% INPUT:
%  Xcal     X-matrix for the calibration set
%  Xval     X-matrix for the validation set
%  Ycal     Y-values for the calibration set
%  Yval     Y-values for the validation set
%  opt      See CV
%
% OUTPUT:
%  model    Predicted Y values

% Copyright, 2014 - 
% This M-file and the code in it belongs to the holder of the
% copyrights and is made public under the following constraints:
% It must not be changed or modified and code cannot be added.
% The file must be regarded as read-only. 
% In case of doubt, contact the holder of the copyrights.
%
% Åsmund Rinnan
% E-mail asmundrinnan@gmail.com

[rY,kY] = size(Ycal);

if length(size(Xcal))~=length(size(Xval)) || size(Xcal,2)~=size(Xval,2) || size(Xcal,3)~=size(Xval,3)
   error('The size of ''Xcal'' and ''Xval'' are not the same')
end
if isempty(Yval)
   Yval = zeros( size( Xval,1), size(Ycal,2));
end
if size(Ycal,2)~=size(Yval,2)
   error('''Ycal'' and ''Yval'' does not have the same number of columns')
end

%Pre-processing of data
%X-matrix
Xcp = Xcal;
Xvp = Xval;
if length( size( Xcal) ) == 2
    for cp = 1:length( opt.xpp.met)
        [Xcp, Xps] = rpls_prepro( Xcp, opt.xpp.met{cp}, opt.xpp.set(cp) );
        Xvp = rpls_prepro( Xvp, opt.xpp.met{cp}, Xps );
    end
    %Remove columns which has been set to NaN (from NW/ SG)
    if any( sum( isnan( Xcp) ) == size( Xcp, 1)) || any( sum( isnan( Xvp) ) == size( Xvp, 1))
        if size( Xvp, 1) == 1
            remid = [sum( isnan( Xcp) ) == size( Xcp, 1); ...
                isnan( Xvp) == 1];
        else
            remid = [sum( isnan( Xcp) ) == size( Xcp, 1); ...
                sum( isnan( Xvp) ) == size( Xvp, 1)];
        end
        remid = sum( remid) > 0;
    else
        remid = false( 1, size( Xcp, 2) );
    end
    
    %Remove columns which do not have any variation in Xcp
    remid( var( Xcp) == 0) = true;
    
    Xcp( :, remid) = [];
    Xvp( :, remid) = [];
    
    model.remid = remid;
else
    [Xcp, Xm, Xs] = nprocess( Xcp, [1 0 0], [0 0 0]); %Only does mean-centering, for now
    Xvp = nprocess( Xvp, [1 0 0], [0 0 0], Xm, Xs);    
end

%Y-matrix
Ycp = Ycal;
Yvp = Yval;
if length( size( Ycal) ) == 2
    for cp = 1:length( opt.ypp.met)
        yv = rpls_nanvar( Ycp);
        model.Yrem = yv == 0 | isnan( yv);
        [Ycp, Yps] = rpls_prepro( Ycp( :, ~model.Yrem), opt.ypp.met{cp}, opt.ypp.set(cp) );
        Yvp = rpls_prepro( Yvp( :, ~model.Yrem), opt.ypp.met{cp}, Yps );
    end
else
    error( 'The code cannot handle multiway data at the moment :(')
end

if isempty( Xcp) %In case all variables has been removed
    model = [];
else    
    switch lower( opt.method)
        case 'pls'
            if sum( rpls_vec( isnan(Xcp) ) ) > 0
                if length( Yps.ref) == 1
                    Yps.ref( 2) = 1;
                end
                for cf = 1:opt.nLV
%                     [b, w, t, p, u, q] = plsmiss( Xcp, Ycp, cf);
                    [mod, xnew] = rpls_plsmiss( Xcp, Ycp, cf);
                    if sum( mod.xid) < size( Xcp, 1)
                        error( 'This function cannot handle samples which are completely missing')
                    end
                    %To find the rms0-value to input into plspred
                    temp = Xcp;
                    tp = mod.Tx * mod.Px';
                    temp( isnan( temp) ) = tp( isnan( temp) );
                    yc = temp * mod.B' * Yps.ref(2) + Yps.ref(1);
                    rms0 = sum( rpls_cen_std( Ycal).^2 );
                    [yp, tp, xr] = rpls_plspred( Xvp, mod.B', mod.Px, Yps.ref, rms0, true); %Not totally flexible, but works for now
                    model.yp( :, cf) = yp( :, end);
                    model.Tp{ cf} = tp;
                    if size( xr, 1) == 1 
                        model.Xres( :, cf) = xr( end); %Changed 131010
                    else
                        model.Xres( :, cf) = xr( :, end);
                    end
                    model.B( :, cf) = mod.B( end, :)';
                    model.Up{ cf} = Yvp * pinv( mod.Py', 3*eps);
                end
            else
%                 [model.B, w, t, p, u, q] = apls( Xcp, Ycp, opt.nLV, [], [], 1);
                if size( Ycp, 2) == 1
                    [model.B, w, t, p, u, q] = rpls_bipls( Xcp, Ycp, opt.nLV);
                else
                    [model.B, w, t, p, u, q] = apls( Xcp, Ycp, opt.nLV);
                end
                [model.yp, model.Tp, model.Xres] = rpls_plspred( Xvp, model.B, p, Yps.ref);
                if sum( rpls_vec( isnan( Xvp) ) ) > 0
                    for cf = 1:length( model.Tp)
                        model.Up{cf} = Yvp * pinv( q( :, 1:cf)', 3 * eps);
                    end
                else
                    model.Up = Yvp * pinv( q', 3*eps);
                end
            end
        case 'pcr'        
        case 'parafac'
        case 'npls'
            [model.Xfac, model.Yfac, model.Core, model.B] = npls( Xcp, Ycp, opt.nLV, NaN);
            for cf = 1:3
                [model.yp( :, cf), model.Tp, ssx, model.Xres{cf}] = npred( Xvp, cf, model.Xfac, model.Yfac, model.Core, model.B, NaN);
            end
        case 'wpls'
            [model.B, w, t, p, u, q] = wpls( Xcp, opt.weights.X, Ycp, opt.weights.Y, opt.nLV);
            [model.yp, model.Tp, model.Xres] = rpls_plspred( Xvp, model.B, p, Yps.ref);
            model.Up = Yvp * pinv( q', 3 * eps);
    end
end
