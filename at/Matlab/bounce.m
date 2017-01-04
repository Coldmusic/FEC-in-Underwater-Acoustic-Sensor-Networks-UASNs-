function bounce( filename )

% run the BOUNCE program
%
% usage: bounce( filename )
% where filename is the environmental file

runbounce = which( 'bounce.exe' );
eval( [ '! "' runbounce '" ' filename ] );


