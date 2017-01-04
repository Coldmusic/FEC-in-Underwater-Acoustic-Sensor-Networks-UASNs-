% Run tests to verify the capability to use a source beam pattern in
% BELLHOP
global units
units = 'm';

bellhop( 'omni' )
plotshd( 'omni.shd', 2, 1, 1 );
caxisrev( [ 50 100 ] )
axis( 'equal' )

bellhop( 'shaded' )
plotshd( 'shaded.shd', 2, 1, 2 );
caxisrev( [ 50 100 ] )
axis( 'equal' )
