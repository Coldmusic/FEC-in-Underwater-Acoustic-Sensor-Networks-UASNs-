function krakenc( filename )

% run the KRAKENC program
%
% usage: krakenc( filename )
% where filename is the environmental file


runkrakenc = which( 'krakenc.exe' );
eval( [ '! "' runkrakenc '" ' filename ] );

runfield = which( 'field.exe' );
eval( [ '! "' runfield '" ' filename ' < field.flp > field.prt' ] );
