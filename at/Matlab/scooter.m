function scooter( filename )

% run the SCOOTER program
%
% usage: scooter( filename )
% where filename is the environmental file

runscooter = which( 'scooter.exe' );
eval( [ '! "' runscooter '" ' filename ] );

% Fortran fields routine:
runfields = which( 'fields.exe' );
eval( [ '! "' runfields '" ' filename ' < fields.flp > fields.prt' ] );

% Matlab fields routine:
% fieldsco( filename );
