function [Yp, Tp, Xres] = rpls_plspred( Xval, B, P, Yscale, rms0, mis)
% Gives you the predicted Y's and T's based on a PLS or PLS-model
%
%[Yp, Tp, Xres] = plspred( Xval, B, P, Yscale, rms0, mis)
%
%INPUT:
% Xval   X used for the prediction (should be pre-processed as the data
%         used in the calibration!)
% B      Regression coefficients for PLS/PCR-model
% P      X-loadings for PLS/PCR-model, if omitted, only Yp in output (only
%         possible for Xval with no NaN-values)
% Ymean  Mean and scaling factors for reference data. < default = [0; 1] >
% rms0   Root-mean-squared error for a 0 factor solution (not necessary, 
%         but good for data with missing values)
% mis    True if the model was based using 'plsmiss', false if not
%
%OUTPUT:
% Yp     Predicted Y
% Tp     Predicted T (score)
% Xres   Residual in X

% 080813 AAR There are cases where the calibration data does NOT contain
%             missing values, but the validation does. Added an extra input
%             value for these cases
% 061210 AAR Improved the prediction when 'Xval' included missing values
% 210908 AAR

if nargin < 3 && sum(vec(isnan(Xval))) > 0
    error('''P'' must be included if your data has missing values')
else
    nLV = size( B, 2);
end

%Sometimes the B came in with the wrong dimensions!
if iscell( B)
    if size( B{1}, 2) == size( Xval, 2)
        for cf = 1:length( B)
            B{cf} = B{cf}';
        end
%         error( 'something is wrong')
    end
end

if ~exist('P', 'var') || isempty( P )
    P = rand( size( Xval, 2), nLV);
end

if ~exist('Yscale', 'var') || isempty( Yscale )
    Yscale = [0; 1];
end

%In case the Y-data only has been mean-centered
if length( Yscale) < 2
    Yscale( 2) = 1;
end

if nargin < 5
    rms0 = 1;
elseif isempty( rms0)
    rms0 = 1;
end

if nargin < 6
    mis = false;
end
% Missing values means that the model are not iterative (?)
nLV = size( P, 2);

e = Xval;
XNaN = isnan( Xval );
if sum(rpls_vec(isnan(Xval))) > 0
    if mis
        st = nLV;
        en = nLV;
    else
        st = 1;
        en = nLV;
    end
    for cf = st:en
        Tp = zeros( size( Xval, 1), cf);
        p = P( :, 1:cf);
        b = B( :, cf);
        Xpred = Xval;
        diff = 10;
        rms = 1;
        e( XNaN) = 0;
        Ypold = zeros( size( Xval, 1), size( b, 2) );
        it = 1;
        while diff > 1e-8 && rms > 1e-4 %Should have q2 rms instead of this number...
            t = e * pinv(p', 3*eps);
            new = t * p';
            diff = sum( (e( XNaN ) - new( XNaN )).^2 );
            e( XNaN ) = new( XNaN );
            Tp = t;
            mod = Tp * p';
            Xpred( XNaN ) = mod( XNaN );
            Yp = (Xpred * b) * Yscale(2) + Yscale(1);
            rms = mean( sqrt( mean( (Ypold - Yp).^2 ) ) )/rms0;
            it = it + 1;
            if it > 1e4
                disp( 'Prediction did not converge')
                break
            end
            Ypold = Yp;
        end
        Ttot{cf} = Tp;
        Yptot( :, cf) = Yp;
    end
    Yp = Yptot;
    if mis
        Tp = Ttot{end};
    else
        Tp = Ttot;
    end
else % No NaN-values
    Tp = Xval * pinv( P', 3*eps);
    if iscell( B)
        for cy = 1:length( B)
            Yp{ cy} = (Xval * B{cy}) * Yscale(2) + Yscale(1);
        end
    else
        Yp = (Xval  * B) * Yscale(2) + Yscale(1);
    end
end

if nargout == 3
    if sum( rpls_vec( isnan( Xval ) ) ) > 0
        if iscell( Tp)
            for cf = 1:length( Tp)
                temp = Tp{cf} * P( :, 1:cf)';
                res = Xval - temp;
                res( XNaN) = 0;
                Xres( :, cf) = sum( res.^2, 2);
            end
        else
            temp = Tp * P';
            res = Xval - temp;
            res( XNaN ) = 0;
            Xres = sum( res.^2, 2);
        end
    else
        Xres = zeros( size( Tp, 1), nLV);
        for cn = 1:nLV
            temp = Tp( :, 1:cn) * P( :, 1:cn)';
            res = Xval - temp;
            Xres(:, cn) = sum( res.^2, 2);
        end
    end
end