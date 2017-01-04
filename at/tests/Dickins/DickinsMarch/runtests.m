% Dickins seamount test
% mbp
global units jkpsflag

units = 'km';

Nprof = 601;  % change field.flp to also
r     = linspace( 0.0, 30000, Nprof );

kraken_rd( 'DickinsK', r )

% Fortran version of field
runfield  = which( 'field.exe' );
eval( [ '! "' runfield '" DickinsK_rd < field.flp' ] );

figure
plotshd( 'DickinsK_rd.shd' )
caxisrev( [ 70 120 ] )

% Matlab version of field
field( 'DickinsK_rd' )

figure
plotshd( 'DickinsK_rd.mat' )
caxisrev( [ 70 120 ] )
