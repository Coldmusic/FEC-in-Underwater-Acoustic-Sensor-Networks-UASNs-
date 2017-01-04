% calibration test cases
% mbp

% isovelocity calibration case

bellhop( 'calibB' )
plotshd( 'calibB.shd', 2, 2, 1 )
caxisrev( [ 40 80 ] )

bellhop( 'calibB_gb' )
plotshd( 'calibB_gb.shd', 2, 2, 2 )
caxisrev( [ 40 80 ] )

kraken( 'calibK' )
plotshd( 'calibK.shd', 2, 2, 3 )
caxisrev( [ 40 80 ] )

scooter( 'calibS' )
plotshd( 'calibS.shd', 2, 2, 4 )
caxisrev( [ 40 80 ] )

% gradient calibration case

bellhop( 'calibBgrad' )
plotshd( 'calibBgrad.shd', 2, 2, 1 )
caxisrev( [ 40 80 ] )

bellhop( 'calibBgrad_gb' )
plotshd( 'calibBgrad_gb.shd', 2, 2, 2 )
caxisrev( [ 40 80 ] )

kraken( 'calibKgrad' )
plotshd( 'calibKgrad.shd', 2, 2, 3 )
caxisrev( [ 40 80 ] )

scooter( 'calibSgrad' )
plotshd( 'calibSgrad.shd', 2, 2, 4 )
caxisrev( [ 40 80 ] )
