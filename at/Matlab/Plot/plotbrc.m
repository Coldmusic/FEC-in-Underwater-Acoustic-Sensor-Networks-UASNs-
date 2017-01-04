function plotbrc( brcfil )

% plot the Bottom Reflection Coefficient file used by Bellhop
% usage: plotbrc( brcfil )
%
% MBP May 2002

global thetaBot RBot phiBot NBotPts thetaTop RTop phiTop NTopPts


% Read the reflection coefficient
trcfil = ' ';
readrc( trcfil, brcfil, ' ', 'F' )

[ AX, ~, H2 ] = plotyy( thetaBot, 20 * log10( RBot ), thetaBot, phiBot, 'plot' );

set( get( AX( 1 ),'Ylabel'), 'String','|R|') 
set( get( AX( 2 ),'Ylabel'), 'String','angle (degrees)' )

xlabel( 'angle (degrees)' )
title( 'Bottom Reflection Coefficient' )

%plot( thetaBot, 20 * log10( RBot ), 'k-' )
%xlabel( 'angle (degrees)' )
%ylabel( 'Reflection Loss (dB)' )
