% Run tests to verify the source beam pattern is incorporated correctly in BELLHOP

global units
units = 'm';

bellhopM( 'omni' )
bellhopM( 'shaded' )

plotshd( 'omni.mat', 2, 1, 1 );
caxisrev( [ 50 100 ] )
axis( 'equal' )

plotshd( 'shaded.mat', 2, 1, 2 );
caxisrev( [ 50 100 ] )
axis( 'equal' )

