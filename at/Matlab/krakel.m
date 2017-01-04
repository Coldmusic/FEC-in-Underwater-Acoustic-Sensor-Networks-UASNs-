function krakel( filename )

% run the KRAKEL program
%
% usage: krakel( filename )
% where filename is the environmental file (without the extension)

runkraken = which( 'krakel.exe' );
eval( [ '! "' runkraken '" ' filename ] );

runfield = which( 'field.exe' );
eval( [ '! "' runfield '" ' filename ' < field.flp > field.prt' ] );

