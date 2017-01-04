function kernelsparc( SSP, c2R, c2I, rho, crosst, cMin, cMax, k, freq, Pos, rr, Bdry, PlotTitle )

% Solve system for a sequence of k-values

global tMult
global Green RTSrd RTSrr tout deltat

% pre-allocate for efficiency
Nrd    = length( Pos.r.depth );
Nr     = length( rr );
Nk     = length( k );
Ntout  = length( tout );
deltak = ( k( end ) - k( 1 ) ) / Nk;

switch( Bdry.Top.Opt(4:4) )
    case ( 'D')
        RTSrd = zeros( Nrd, Ntout );
    case ( 'R')
        RTSrr = zeros( Nr,  Ntout );
    otherwise
        Green = zeros( Ntout, Nrd, Nk );
end

deltat = tMult / sqrt( 1.0 / crosst^2 + ( 0.5 * cMax * k( Nk ) )^2 ); % Courant condition to set time step

fprintf( '\n\nTime step (based on CFL condition) = %f', deltat )
fprintf( '\n\nEstimated fl. pt. ops (millions) = %f \n', ( tout( end ) / deltat ) * Nk * length( c2R ) / 25000 )

% *** Loop over spectral components ***

for Ik = 1 : Nk
    x = k( Ik )^2;
    deltat = tMult / sqrt( 1.0 / crosst^2 + ( 0.5 * cMax * k( Ik ) )^2 );
    if ( mod( Ik, 10 ) == 0 )   % print every 10th wavenumber
       fprintf( '\n Ik, Nk %i %i', Ik, Nk )
    end
    march( SSP, c2R, c2I, rho, x, Ik, deltak, rr );  % March that component for all time
end

switch ( Bdry.Top.Opt(4:4) )
    case ( 'S' )          % *** Case of a snapshot ***
        PlotType = 'rectilin  ';
        Pos.r.range = k;    % store wavenumbers in slot usually used for receiver ranges
        Pos.s.depth = tout;  % store times       in slot usual used for source depths
        p(1, :, :, : ) = Green;  % store Green's function in slot usually used for pressure
        atten = 0;
        save GRNFIL PlotTitle PlotType freq atten Pos p

    case ( 'D' )   % *** Write out RTS (vertical array) *
        Scale = 1.0 / sqrt( pi * rr( 1 ) );
        Pos.r.range = k;   % store wavenumbers is slot usually used for receiver ranges
        RTS = Scale * RTSrd;
        save RTSFIL PlotTitle Pos tout RTS

    case ( 'R' )   %  *** Write out RTS (horizontal array) ***
        Pos.r.range = k;   % store wavenumbers is slot usually used for receiver ranges
        RTS = RTSrd;
        save RTSFIL PlotTitle Pos tout RTS
end
