function [ pltitl, plottype, freq, atten, Pos, p ] = read_shd_asc( filename )

% Read an ascii shade file

% open the file
fid = fopen( filename, 'r' );
if ( fid == -1 )
    errordlg( 'No shade file with that name exists; you must run a model first', 'read_shd_asc' );
    error(    'No shade file with that name exists; you must run a model first', 'read_shd_asc' );
end

% read

pltitl   = fgetl( fid );
plottype = fgetl( fid );
freq     = fscanf( fid, '%f', 1 );
Ntheta   = fscanf( fid, '%i', 1 );
Nsd      = fscanf( fid, '%i', 1 );
Nrd      = fscanf( fid, '%i', 1 );
Nrr      = fscanf( fid, '%i', 1 );
atten    = fscanf( fid, '%f', 1 );

Pos.theta    = fscanf( fid, '%f', Ntheta );
Pos.s.depth  = fscanf( fid, '%f', Nsd );
Pos.r.depth  = fscanf( fid, '%f', Nrd );
Pos.r.range  = fscanf( fid, '%f', Nrr );

isd = 1;
for i = 1:isd
   temp1   = fscanf( fid, '%f', [ 2 * Nrr, Nrd ] );
   i, size(temp1)
end

fclose( fid );

% joint real and imaginary parts into a complex matrix

p = temp1( 1:2:2*Nrr, : )' + sqrt( -1 ) * temp1( 2:2:2*Nrr, : )';