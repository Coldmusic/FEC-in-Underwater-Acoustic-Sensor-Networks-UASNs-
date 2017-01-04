function reflect( I, BeamType, BC, cHS, rhoHS, BotTop, TBdry, NBdry, kappa, SSP )
% Reflect a ray/beam off a boundary

global ray omega Layer

RADDEG = 180 / pi;   % used to convert radians to degrees

Tg = dot( ray( I ).Tray', TBdry' );  % component of ray tangent, along boundary
Th = dot( ray( I ).Tray', NBdry' );  % component of ray tangent, normal to boundary

% *** calculate the change in curvature ***
% Based on formulas given by Muller, Geoph. J. R.A.S., 79 (1984).

[ c, gradc, ~, ~, ~, Layer ] = ssp( ray( I ).x( : ), SSP, Layer );

cn = gradc( 2 ) .* ray( I ).Tray( 1 );
cs = gradc( 2 ) .* ray( I ).Tray( 2 );   % assumes gradc( 2 ) = cr = 0

RN = 2 * kappa / c^2 / Th;    % boundary curvature correction

if ( strcmp( BotTop, 'TOP' ) )
    cn = -cn;    % flip sign for top reflection
    RN = -RN;
end

RM = Tg ./ Th;
RN = RN + RM' .* ( 4 * cn - 2 * RM' .* cs ) ./ c;


switch ( BeamType(2:2) )
    case ( 'D' )
        RN = 2.0 * RN;
    case ( 'Z' )
        RN = 0.0;
end

ray( I ).Tray( 1 ) =  ray( I ).Tray( 1 ) - 2.0 * ( Th .* NBdry( 1 )' )';
ray( I ).Tray( 2 ) =  ray( I ).Tray( 2 ) - 2.0 * ( Th .* NBdry( 2 )' )';

ray( I ).p( 1 ) = ray( I ).p( 1 ) + ray( I ).q( 1 ) * RN;
ray( I ).p( 2 ) = ray( I ).p( 2 ) + ray( I ).q( 2 ) * RN;

% *** account for phase change ***

switch ( BC )
    case ( 'R' )                 % rigid
    case ( 'V' )                 % vacuum
        ray( I ).Rfa = -ray( I ).Rfa;
    case ( 'F' )                 % file
        theInt = RADDEG * abs( atan2( Th, Tg ) );   % angle of incidence (relative to normal to bathymetry)
        %[ theta, RefC, phi ] = RefCO( theInt, rInt, phiInt, Npts  )
        ray( I ).Rfa = ray( I ).Rfa * rInt * exp( 1i * phiInt );
    case ( 'A' )                 % half-space
        GK       = omega * Tg;   % wavenumber in direction parallel to bathymetry
        gamma1SQ = ( omega / c'  ).^ 2 - GK^ 2;
        gamma2SQ = ( omega / cHS ) ^ 2 - GK^ 2;
        gamma1   = sqrt( -gamma1SQ );
        gamma2   = sqrt( -gamma2SQ );
        Refl = ( rhoHS * gamma1 - gamma2 ) / ( rhoHS * gamma1 + gamma2 );
        ray( I ).Rfa = Refl.' .* ray( I ).Rfa;

end
