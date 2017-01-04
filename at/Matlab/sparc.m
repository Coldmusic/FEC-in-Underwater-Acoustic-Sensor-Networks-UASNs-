function sparc( filename )

% run the SPARC program
%
% usage: sparc( filename )
% where filename is the environmental file

runsparc = which( 'sparc.exe' );
eval( [ '! "' runsparc '" ' filename ] );
