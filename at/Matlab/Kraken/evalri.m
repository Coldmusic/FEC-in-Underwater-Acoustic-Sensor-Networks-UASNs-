function p = evalri( phiS, phi, R, k, Opt, Comp )

% conversion from Fortran eval.f90
% mbp 4/2009
%
%
% Computes pressure field using coupled mode theory
% Normalized to pressure of point source at 1 meter
%
% Opt = X     Cartesian   (X, z) coordinates
% Opt = R     Cylindrical (R, z) coordinates

% For vertical component of displacement field, take finite difference in depth
% This is done using a backward difference.
% A centered difference or an FFT formula would be better

if ( Comp == 'V' )
   phidiff = diff( phi );   % should divide by diff( z ) also
   phi( 1: end - 1, : ) = phidiff;
   phi( end, : )        = zeros( 1, size( phi, 2 ) );   % no derivative for last row
end

phi = phi * diag( phiS, 0 );	% scale modes by phiS

% form pressure field

phase = diag( 1.0 ./ sqrt( k ) ) * exp( -1i * k * R' ) * diag( realsqrt( 2 * pi ./ R ) );

% for horizontal component take derivative in range direction
% The following formula approximates the derivate, assuming the phase term
% (e^(-i k rr )) dominates
if ( Comp == 'H' )
   phase = diag( -1i .* sqrt( k ) ) * exp( -1i * k * R' ) * diag( realsqrt( 2 * pi ./ R ) );
end

p( 1, 1, :, : ) = phi * phase;
