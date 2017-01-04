function kraken( filename )

% run the KRAKEN program
%
% usage: kraken( filename )
% where filename is the environmental file (without the extension)

runkraken = which( 'kraken.exe' );
eval( [ '! "' runkraken '" ' filename ] );

runfield = which( 'field.exe' );
eval( [ '! "' runfield '" ' filename ' < field.flp > field.prt' ] );

