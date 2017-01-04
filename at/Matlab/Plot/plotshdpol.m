function plotshdpol( filename )

% plot a single TL surface in dB (polar coordinates)
% usage: plotshdpol( filename )

% open the file and read data

[ PlotTitle, ~, ~, ~, Pos, pressure ] = read_shd( filename );
pressure = squeeze( pressure( :, 1, 1, : ) );   % take first source and receiver depth

theta = Pos.theta;
rkm   = Pos.r.range / 1000.0;

% shift coordinate system to arbitrary position
xs = 333;	% -12.0
ys = 315;	%  13.8

% make plot polar

[th, r ] = meshgrid( theta, rkm );
th = ( 2 * pi / 360. ) * th;   % convert to radians
[x, y] = pol2cart( th, r );

x = x + xs * ones( size( x ) );
y = y + ys * ones( size( x ) );

tlt = -20.0 * log10( abs( pressure ) );
tlt = tlt';

% *** plot ***

tej = flipud( colormap( jet( 256 ) ) );
figure
surfc( x, y, tlt ); shading interp; ...
colormap( tej ); caxis( [ 75 95 ] ); colorbar; view( 2 )
xlabel( 'Range (km)' ); ylabel( 'Range (km)' );
title( deblank( PlotTitle ) )
axis( 'equal' )