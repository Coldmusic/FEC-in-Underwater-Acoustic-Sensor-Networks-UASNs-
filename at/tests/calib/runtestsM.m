% calibration test cases
% mbp

% isovelocity calibration case

bellhopM( 'calibB' )
plotshd( 'calibB.mat', 2, 2, 1 )
caxisrev( [ 40 80 ] );

bellhopM( 'calibB_gb' )
plotshd( 'calibB_gb.mat', 2, 2, 2 )
caxisrev( [ 40 80 ] );

kraken( 'calibK' )
plotshd( 'calibK.shd', 2, 2, 3 )
caxisrev( [ 40 80 ] );

scooterM( 'calibS' )
plotshd( 'calibS.mat', 2, 2, 4 )
caxisrev( [ 40 80 ] );

% gradient calibration case

bellhopM( 'calibBgrad' )
plotshd( 'calibBgrad.mat', 2, 2, 1 )
caxisrev( [ 40 80 ] );

bellhopM( 'calibBgrad_gb' )
plotshd( 'calibBgrad_gb.mat', 2, 2, 2 )
caxisrev( [ 40 80 ] );

kraken( 'calibKgrad' )
plotshd( 'calibKgrad.shd', 2, 2, 3 )
caxisrev( [ 40 80 ] );

scooterM( 'calibSgrad' )
plotshd( 'calibSgrad.mat', 2, 2, 4 )
caxisrev( [ 40 80 ] );
