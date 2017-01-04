function [ Modes ] = read_modes_asc( filename, modes )

% Read the modes produced by KRAKEN
% usage:
%    [ pltitl, freq, ck, z, phi ] = read_modes_asc( filename, modes )
% filename should include extension
% modes is an optional vector of mode indices

fid = fopen( filename, 'r' );

if ( fid == -1 )
   errordlg( 'Mode file does not exist', 'read_modes_asc' )
   error(    'read_modes_asc: Mode file does not exist' )
end

lrecl   = fscanf( fid, '%i' );
Modes.pltitl = fgetl( fid );

temp = fscanf( fid, '%f', 5 );
Modes.freq   = temp( 1 );
Modes.Nmedia = temp( 2 );
Modes.ntot   = temp( 3 );
Modes.nmat   = temp( 4 );
Modes.M      = temp( 5 );

%temp   = fscanf( fid, '%i', [ nmedia ] )
for ii = 1 : Modes.Nmedia
   junk   = fgetl( fid );
end
junk   = fgetl( fid );   % top halfspace properties
junk   = fgetl( fid );   % bot halfspace properties
junk   = fgetl( fid );   % blank line

Modes.z      = fscanf( fid, '%f', [ 1, ntot ] );
junk   = fgetl( fid );
ckt    = fscanf( fid, '%f', [ 2, m ] );
Modes.ck = ckt( 1, : )' + 1i * ckt( 2, : )';

if nargin == 1
   modes = 1:m;    % read all modes if the user didn't specify
end

% don't try to read modes that don't exist
ii =  modes <= m ;
modes = modes( ii );

Modes.phi = zeros( ntot, length( modes ) );

for mode = 1: max( modes )
   junk = fgetl( fid );
   phit = fscanf( fid, '%f', [ 2, ntot] );  % read the mode
   
   ii = find( mode == modes );  % see if it's in the list to grab
   if ii >= 1                   % if yes, the store it
      Modes.phi( :, ii )  = phit( 1, : )' + 1i * phit( 2, : )';
   end
end

%figure; imagesc( 1:m, z, real( phi ) ); colorbar
%figure; surf( 1:m, z, real( phi ) ); colorbar

