function model = rpls_simpls( X, Y, nLV)
%model = rpls_simpls( X, Y, nLV)
%
%INPUT:
% X      Matrix of independent variables
% Y      Matrix of dependent variables
% nLV    Number of factors to be extracted
%
%OUTPUT:
% model  A structure with model parameters

% AAR 051118 Please note that the NIPALS and SIMPLS does NOT give the same
%             result for the PLS2 case. They're similar, but not identical
% AAR 030810 Based on Jong S., Chemometrics and Intelligent Laboratory 
%             Systems, 18, 251-263 and Andersson M., Journal of Chemometrics, 
%             23, 518-529

covXY = X' * Y;
[rX, cX] = size( X);
[rY, cY] = size( Y);

if rX ~= rY
    error( 'The number of samples in X and Y must be the same')
end

%Initializing the matrices
Tx = zeros( rX, nLV);
Px = zeros( cX, nLV);
Ty = zeros( rY, nLV);
Py = zeros( cY, nLV);
basis = Px;

for cf = 1:nLV
    
    %The SVD of the residual covariance matrix gives the y-loading
    [u, s, v] = svds( covXY' * covXY, 1);
    pY = u;
    %Multiplying this with covXY and normalizing gives the loading weights
    w = covXY * pY;        
    w = normalize( w); %091116 New line
    
    %Calculate the X-scores from X and the calculated loading weights
    tx = X * w; %tt
    
    %091116 These lines have been deleted
%         nt = norm( tx); %normt 
%     tx = tx/ nt; %tt
%     w = w * nt; %rr X = tx * px' -> X' = px * tx' -> X' * tx = px * tx' *
%     tx -> px = X' * tx * inv( tx' * tx);

    %Calculate X-loadings from X and the X-scores
    px = X' * tx * inv( tx' * tx); %(tx' * X)'; %pp

%     %091116 New line
%     px = X' * tx * inv( tx' * tx);    
%     
    %Calculate the Y-loadings from Y and the X-scores
    py = Y' * tx * inv( tx' * tx); %qq   Y = tx * py' -> py = Y' * tx * inv( tx' * tx)
    
%     %091116 New line
%     py = normalize( py);
    
    %Calculate the Y-scores from the Y and the Y-loadings.
    ty = Y * py * inv( py' * py); %uu
    
    %Set v equal to the alst Y-loading
    v = px; %vv
    if cf > 1
        %Remove anything in 'v' that has been explained earlier. ('basis'
        %is the rest of the model
        v = v - basis * (basis' * px);
        %Do the same for the Y-score where the X-scores are the rest of the
        %model
        ty = ty - Tx * (Tx' * ty);
    end
    v = normalize( v);
    covXY = covXY - v * inv( v' * v) * (v' * covXY);
    %Reduce the covariance matrix accordingly. 
%     covXY = covXY - v * (v' * covXY);
    W( :, cf) = w;
    Tx( :, cf) = tx;
    Px( :, cf) = px;
    Ty( :, cf) = ty;
    Py( :, cf) = py;
    basis( :, cf) = v;
end
model.Tx = Tx;
model.Px = Px;
model.Ty = Ty;
model.Py = Py;
model.W = W;
model.Xbase = basis;
        
%Let's normalize so that W is orthoNORMAL and the length goes into Tx

model.B = cell( cY, 1);
for ny = 1:cY
    if size( W, 2) > 1
        model.B{ny} = cumsum( ( W * diag( Py( ny, :) ) )' )';
    else
        model.B{ny} = ( ( W * Py( ny) )')';
    end
end

if cY == 1
    model.B = model.B{1};
end