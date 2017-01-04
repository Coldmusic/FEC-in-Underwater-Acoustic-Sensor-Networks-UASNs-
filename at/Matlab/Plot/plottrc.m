function plottrc( trcfil )

% plot the Top Reflection Coefficient file used by Bellhop
% usage: plottrc( trcfil )
%
% MBP May 2002

global thetaBot RBot phiBot NBotPts thetaTop RTop phiTop NTopPts

% Read the reflection coefficient
brcfil = ' ';
readrc( trcfil, brcfil, 'F', ' ' )

[ AX, ~, H2 ] = plotyy( thetaTop, RTop, thetaTop, phiTop, 'plot' );

set( get( AX( 1 ),'Ylabel'), 'String','|R|') 
set( get( AX( 2 ),'Ylabel'), 'String','angle (degrees)' )

xlabel( 'angle (degrees)' )
title( 'Top Reflection Coefficient' )

%plot( thetaBot, 20 * log10( RBot ), 'k-' )
%xlabel( 'angle (degrees)' )
%ylabel( 'Reflection Loss (dB)' )
