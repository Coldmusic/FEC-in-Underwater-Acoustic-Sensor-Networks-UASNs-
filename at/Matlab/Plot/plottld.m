function plottld( filename, rrt )

% plot a single TL slice from the shade file
%
% usage:
% plottld( filename, rrt )
% where
%   filename is the shadefile (with extension)
%   rrt is the receiver range in km
%   if rrt is a vector then one plot is generated for each element
% mbp

% read

[ PlotTitle, ~, ~, ~, Pos, pressure ] = read_shd( filename );

disp( 'PlotTLd uses the first bearing and source depth in the shade file; check OK' )
itheta = 1
isd    = 1

tlt = abs( pressure );	            % this is really the negative of TL
%ii = find( tlt == 0 );                      % find places where pressure vanishes
%tlt( ii ) = max( max( tlt ) ) / 1e10;       % and set it to a small number
tlt( tlt == 0 ) = max( max( tlt ) ) / 1e10;  % find places where pressure vanishes and set it to a small number

tlt = -20.0 * log10( tlt );          % so there's no error when we take the log

% following logic is because interp1 won't interpolate a vector with only 1 element
if ( length( Pos.r.depth ) == 1 )
  tlslice = tlt;
else
  TLtemp = squeeze( tlt( itheta, isd, :, : ) )';   % need to avoid dimensional problems when TLT has 1, 2, or 3 dimensions
  if ( length( Pos.r.range ) == 1 )
     tlslice = TLtemp( : )';
  else
     tlslice = interp1( Pos.r.range, TLtemp, 1000.0 * rrt );
  end
end

hh= plot( tlslice, Pos.r.depth )
set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
set( gca, 'Xdir', 'Reverse' )
xlabel( 'TL (dB)'   )
ylabel( 'Depth (m)' )
title( deblank( PlotTitle ) )
set( hh, 'LineWidth', 2 )

% generate legend
for irr = 1: length( rrt )
    legendstr( irr, : ) = [ 'Range = ', num2str( rrt( irr ) ), ' km' ];
end

legend( legendstr, 'Location', 'Best' )
legend( 'boxoff' )
drawnow
