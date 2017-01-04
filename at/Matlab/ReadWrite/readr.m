function [ rr ] = readr( fid )

% Read receiver ranges

Nrr = fscanf( fid, '%i', 1 );
fprintf( '\nNumber of receiver ranges = %i \n', Nrr )
fgetl( fid );

rr = fscanf( fid, '%f', Nrr );

fprintf( '\nReceiver ranges (m) \n' )
fprintf( '%f ', rr )
fprintf( '\n' )

if Nrr > 2
  rr = linspace( rr( 1 ), rr( 2 ), Nrr )'; % generate vector of receiver ranges
end

fgetl( fid );
