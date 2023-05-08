function [xpp, ppset] = rpls_prepro( x, met, opt)
%[xpp, ppset] = rpls_prepro( x, met, opt)
%
% Pre-processing of spectra. So far the same pre-processing is used on the
% whole spectra
%
%INPUT:
% x     Spectral data with samples as rows
% met   Method - SNV, MSC, SG (Savitzky-Golay), MC (mean-centering), AS
%        (autoscaling) and PS (Pareto-scaling)
% opt   Settings for the pre-processing wanted. Number of inputs depends on
%        method
%
%OUTPUT:
% xpp   Pre-processed spectra
% ppset Values needed for pre-processing, i.e. mean( xcal) for MSC
%
%See also: MSCT, SNV, SG

% 070110 AAR

if nargin == 0
    xpp = struct( 'ref', NaN, 'ax', NaN, 'set', NaN);
    return
end

if ~isstruct( opt)
    opt = rpls_prepro;
end

if nargin < 3
    error( '''opt'' must be given even if it is empty (as for SNV)')
end

ppset = opt;

met = lower( met);
if sum( strcmp( met, {'snv', 'msc', 'sg', 'mc', 'as', 'ps'}) ) == 0
    error( 'Method not recognized. Choose between snv, msc, sg, mc, as and ps')
end

if iscell( met)
    met = met{1};
end

switch met
    case 'snv'
        xpp = snv( x); %Could be enhanced with including MAD/ Quartile range, etc
    case 'msc'
        if isnan( opt.ref)
            opt.ref = mean( x);
        end
        if isnan( opt.ax)
            ax = 1:size( x, 2);
        end
        if isnan( opt.set)
            opt.set = [1 2]; %EMSC
        end
        xpp = rpls_msct( x, opt.ref, opt.ax, opt.set(1), opt.set(2), NaN, opt.set(3));
        ppset.ref = opt.ref;
        ppset.ax = opt.ax;
    case 'sg'
        xpp = rpls_sg( x, opt.set(1), opt.set(2), opt.set(3), 0);
    case 'mc'
        if isnan( opt.ref)
            [xpp, ppset.ref] = rpls_cen_std( x);
        else
            xpp = rpls_cen_std( x, opt.ref);
        end
    case 'as'
        if isnan( opt.ref)
            [xpp, xm, xs] = rpls_cen_std( x);
            ppset.ref = [xm; xs];
        else
            xpp = rpls_cen_std( x, opt.ref( 1, :), opt.ref( 2, :) );
        end
    case 'ps'
        if isnan( opt.ref)
            ppset.ref = [rpls_nanmean( x); sqrt( rpls_nanstd( x) )];
            xpp = rpls_cen_std( x, ppset.ref( 1, :), ppset.ref( 2, :) );
        else
            xpp = rpls_cen_std( x, opt.ref( 1, :), opt.ref( 2, :) );
        end
end