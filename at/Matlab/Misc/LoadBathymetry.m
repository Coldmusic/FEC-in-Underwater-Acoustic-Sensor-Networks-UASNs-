function [ Bathy ] = LoadBathymetry( xyzfile )
% Loads the bathymetry file
% assumed to be a rectangular grid of data written in the standard GEODAS
% xyz form
%
% Returns a structure:
% Bathy.depth depths
% Bathy.Lat   latitudes
% Bathy.Lon   longitudes
% Bathy.Latkm mapping of lat to km assuming (x,y)=(0,0) corresponds to ( Lon(1), Lat(1) )
% Bathy.Lonkm
%
% mike porter 2/2009 based on earlier version by D. Hernandez

Bathy.fname = xyzfile;
xyz         = load( xyzfile );

Lat  = unique( xyz( :, 2 ) ); % get unique latitudes
Lon  = unique( xyz( :, 1 ) ); % get unique longitudes
Lat  = fliplr( Lat );

nlon = length( Lon );
nlat = length( Lat );

z = reshape( xyz( :, 3 ), nlon, nlat );
z = z';
z = flipud( z );

z( z > 0 ) = NaN;   % remove land

Bathy.depth = -z;
Bathy.Lat   = Lat;
Bathy.Lon   = Lon;

% Create grid in km
[ LatTotkm, ~, A21 ] = dist_wh( [Lat( 1 ) Lat( end ) ] , [ Lon( 1 ) Lon( 1   ) ] );
[ LonTotkm, ~, A21 ] = dist_wh( [Lat( 1 ) Lat( 1   ) ] , [ Lon( 1 ) Lon( end ) ] );

% mapping of lat/long to km
Bathy.Latkm = linspace( 0, LatTotkm, length( Lat ) ) / 1000;
Bathy.Lonkm = linspace( 0, LonTotkm, length( Lon ) ) / 1000;
