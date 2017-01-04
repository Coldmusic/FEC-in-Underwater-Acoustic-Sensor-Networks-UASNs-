function field( FileRoot, modes )

% calculates the field using modes produced by KRAKEN
% parameters read from field.flp
%
% usage: field( filename, modes )
% mbp

% Comp     = 'P';   % select component (P, H, V) (Pressure, Horizontal, Vertical)

clear read_modes_bin % to force rewind to beginning of mode file

[ TitleEnv, Opt, Comp, MLimit, NProf, rProf, R, Pos ] = read_flp;

rd  = Pos.r.depth;   % receiver depths
Nrd = length( rd );
Nr  = length( R );   % receiver ranges

if ( NProf == 1 )               % Range-independent case
   
   filename = [ FileRoot '.mod' ];
   
   if nargin == 1
      Modes = read_modes( filename );
   else
      Modes = read_modes( filename, modes );
   end
   MSrc = length( Modes.k );
   M    = min( MLimit, MSrc );        % Set number of propagating modes
   
   % weight modes by mode excitation
   zs  = Pos.s.depth;
   isd = find( Modes.z >= zs );    % index of source depth
   isd = isd( 1 );
   C   = Modes.phi( isd, : ).';    % column vector with modal weights
   
   % calculate mode values at receiver depths
   irdvec = zeros( 1, length( Pos.r.depth ) );
   for ii = 1 : length( Pos.r.depth )
      zr  = Pos.r.depth( ii );
      ird = find( Modes.z >= zr );	% index of source depth
      irdvec( ii ) = ird( 1 );
   end
   phiR = Modes.phi( irdvec, : );
   
   p = evalri( C, phiR, R, Modes.k, Opt, Comp );
else
   if ( Opt(2:2) == 'C' )       % Range-dependent case
      % Coupled mode theory
      clear evalcm % to force re-initialization of profile ranges
      p = evalcm( FileRoot, Pos, rProf, NProf, rd, Nrd, R, Nr, MLimit, Opt );
   else
      % Adiabatic mode theory
      p = evalad( FileRoot, Pos, rProf, NProf, rd, Nrd, R, Nr, M, Opt );
   end
end

PlotTitle   = TitleEnv;
PlotType    = 'rectilin  ';
freq        = 0.0; % Modes.freq;
atten       = 0.0;
Pos.r.range = R;

save( FileRoot, 'PlotTitle', 'PlotType', 'freq', 'atten', 'Pos', 'p' )
