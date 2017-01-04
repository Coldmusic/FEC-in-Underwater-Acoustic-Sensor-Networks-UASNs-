function RAMtoSHD
%
filename = 'tl.grid';
fid = fopen( filename, 'rb' );

dz    = 0.2;
zmplt = 100;
dr    = 9.98;
rmax  = 100000;

dz = 2.99;
zmplt = 1500;

% read in the grid size information from the header
dz    = fread( fid, 1, 'float32' );   % not sure why this extra word has to be read
dz    = fread( fid, 1, 'float32' );
ndz   = fread( fid, 1, 'int32'   );
zmplt = fread( fid, 1, 'float32' );
dr    = fread( fid, 1, 'float32' );
ndr   = fread( fid, 1, 'int32'   );
rmax  = fread( fid, 1, 'float32' );
freq  = fread( fid, 1, 'float32' );
zs    = fread( fid, 1, 'float32' );
fseek( fid, 56, 'bof' ); %reposition to end of first record

Pos.s.depth = zs;
Pos.r.depth = 0 : ndz * dz : zmplt; % + ndz * dz;
Pos.r.range = 0 : ndr * dr : rmax;

Nrd = length( Pos.r.depth );
Nrr = length( Pos.r.range );

if ( fid == -1 )
   error( 'No shade file with that name exists; you must run a model first' );
end

TL = fread( fid, [ Nrd, Nrr ], 'float32' );    %Read complex data

fclose( fid );

PlotTitle = 'RAM';
PlotType  = 'rectilin  ';
atten     = 0;
p( 1, 1, :, : ) = 10 .^ ( -TL / 20 );
save SHDFIL PlotTitle PlotType freq atten Pos p  % problem: sd is outside the loop
