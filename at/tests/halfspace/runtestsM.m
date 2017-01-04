% Run tests to verify the surface and bottom reflection coefficents are incorporated correctly

scooterM( 'lower_halfS' )
plotshd( 'lower_halfS.mat', 2, 2, 1 );
caxisrev( [ 40 80 ] )

bellhopM( 'lower_halfB' )
plotshd( 'lower_halfB.mat', 2, 2, 2 );
caxisrev( [ 40 80 ] )

scooterM( 'upper_halfS' )
plotshd( 'upper_halfS.mat', 2, 2, 3 );
caxisrev( [ 40 80 ] )

bellhopM( 'upper_halfB' )
plotshd( 'upper_halfB.mat', 2, 2, 4 );
caxisrev( [ 40 80 ] )
