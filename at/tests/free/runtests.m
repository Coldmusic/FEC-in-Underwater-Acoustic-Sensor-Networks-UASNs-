% Run tests to verify the point vs. line source options

% BELLHOP point source cases
bellhop 'freePointB'
plotshd( 'freePointB.shd', 3, 2, 1 );
caxisrev( [ 40 80 ] )

bellhop 'freePoint_gbtB'
plotshd( 'freePoint_gbtB.shd', 3, 2, 3 );
caxisrev( [ 40 80 ] )

bellhop 'freePoint_CervenyB'
plotshd( 'freePoint_CervenyB.shd', 3, 2, 5 );
caxisrev( [ 40 80 ] )

%%

% BELLHOP line source cases
bellhop 'freeLineB'
plotshd( 'freeLineB.shd', 3, 2, 2 );
caxisrev( [ 0 25 ] )

bellhop 'freeLine_gbtB'
plotshd( 'freeLine_gbtB.shd', 3, 2, 4 );
caxisrev( [ 0 25 ] )

bellhop 'freeLine_CervenyB'
plotshd( 'freeLine_CervenyB.shd', 3, 2, 6 );
caxisrev( [ 0 25 ] )

%%

% SCOOTER

copyfile( 'fieldsPoint.flp', 'fields.flp' );
scooter 'freeS'
plotshd( 'freeS.shd', 2, 2, 1 );
caxisrev( [ 40 80 ] )

copyfile( 'fieldsLine.flp', 'fields.flp' )
scooter 'freeS'
plotshd( 'freeS.shd', 2, 2, 2 );
caxisrev( [ 0 25 ] )

%%
% exact (analytic) solution

freq = 5;
k0 = 2 * pi * freq / 1500.0;
zs = 3000;
zr = 0:25:5000;
rr = 25:25:10000;
for iz = 1: length( zr )
   rmat1 = sqrt( rr.^2 + (zr( iz ) - zs )^2 );
   rmat2 = sqrt( rr.^2 + (zr( iz ) + zs )^2 );
   ppoint( iz, : ) = exp( i * k0 * rmat1 ) ./ rmat1;
   pline( iz, : ) = besselh( 0, k0 * rmat1 );
end

pline = sqrt( pi / 2 ) * pline; %normalization to roughly match acoustics toolbox

subplot( 2, 2, 3 )
pcolor( rr, zr, -20 * log10( abs( ppoint ) ) ); shading interp; view( 0, -90 )
caxisrev( [ 40 80 ] )
title( 'Exact point source solution' )

subplot( 2, 2, 4 )
pcolor( rr, zr, -20 * log10( abs( pline ) ) ); shading interp; view( 0, -90 )
caxisrev( [  0 25 ] )
title( 'Exact line source solution' )

