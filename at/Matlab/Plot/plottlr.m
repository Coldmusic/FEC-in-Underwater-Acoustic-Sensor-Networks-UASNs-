function plottlr( filename, rdt )

% plot a single TL slice from the shade file
%
% usage:
% plottlr( filename, rdt )
% where
%   filename is the shadefile (with extension)
%   rdt is the receiver depth in m
%   if rdt is a vector then one plot is generated for each element
% mbp

global units

% read

disp( 'PlotTLr uses the first bearing and source depth in the shade file; check OK' )
itheta = 1
isd    = 1

[ PlotTitle, ~, ~, ~, Pos, pressure ] = read_shd( filename );

rkm = Pos.r.range / 1000.0;          % convert to km
tlt = abs( pressure );	            % this is really the negative of TL
tlt( tlt == 0 ) = max( max( tlt ) ) / 1e10;      % replaces zero by a small number
tlt = -20.0 * log10( tlt );          % so there's no error when we take the log

% following logic is because interp1 won't interpolate a vector with only 1 element

if ( length( Pos.r.depth ) == 1 )
  tlslice = squeeze( tlt( itheta, isd, 1, : ) );
else
  tlslice = interp1( Pos.r.depth, squeeze( tlt( itheta, isd, :, : ) ), rdt );
end

hh = plot( rkm, tlslice' );
set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
xlabel( 'Range (km)' )
ylabel( 'TL (dB)' )
title( deblank( PlotTitle ) )
set( hh, 'LineWidth', 2 )

% generate legend
for ird = 1: length( rdt )
    legendstr( ird, : ) = [ 'Depth = ', num2str( rdt( ird ) ), ' m' ];
end

legend( legendstr, 'Location', 'Best' )
legend( 'boxoff' )
drawnow
