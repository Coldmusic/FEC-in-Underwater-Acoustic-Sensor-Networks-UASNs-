function [ title, plottype, freq, atten, Pos, pressure ] = read_shd_bin( filename )

% Read TL surfaces from a binary Bellhop/Kraken .SHD file
% without having to convert to ASCII first.
% Output is a 4-D pressure field p( Ntheta, Nsd, Nrd, Nrr )
% Chris Tiemann, Feb. 2001

fid = fopen( filename, 'rb' );
if ( fid == -1 )
    errordlg( 'No shade file with that name exists; you must run a model first', 'read_shd_bin' );
    error(    'read_shd_bin: No shade file with that name exists; you must run a model first' );
end

recl     = fread( fid,  1, 'int32' );     %record length in bytes will be 4*recl
title    = fread( fid, 80, '*char' )';

fseek(fid, 4*recl, -1); %reposition to end of first record
plottype = fread( fid, 10, '*char'   );
xs       = fread( fid,  1, 'float32' );
ys       = fread( fid,  1, 'float32' );
plottype = plottype';

fseek(fid, 2 * 4 * recl, -1 ); %reposition to end of second record
freq   = fread( fid, 1, 'float32' );
Ntheta = fread( fid, 1, 'int32'   );
Nsd    = fread( fid, 1, 'int32'   );
Nrd    = fread( fid, 1, 'int32'   );
Nrr    = fread( fid, 1, 'int32'   );
atten  = fread( fid, 1, 'float32' );

fseek( fid, 3 * 4 * recl, -1 ); %reposition to end of record 3
Pos.theta   = fread( fid, Ntheta, 'float32' );

fseek( fid, 4 * 4 * recl, -1 ); %reposition to end of record 4
Pos.s.depth = fread( fid, Nsd,    'float32' );

fseek( fid, 5 * 4 * recl, -1 ); %reposition to end of record 5
Pos.r.depth = fread( fid, Nrd,    'float32' );

fseek( fid, 6 * 4 * recl, -1 ); %reposition to end of record 6
Pos.r.range = fread( fid, Nrr,    'float32' );

%Each record holds data from one source depth/receiver depth pair

switch plottype
    case 'rectilin  '
        pressure = zeros( Ntheta, Nsd, Nrd, Nrr );
        Nrcvrs_per_range = Nrd;
    case 'irregular '
        pressure = zeros( Ntheta, Nsd, 1, Nrr );
        Nrcvrs_per_range = 1;
    otherwise
        pressure = zeros( Ntheta, Nsd, Nrd, Nrr );
        Nrcvrs_per_range = Nrd;
end

for itheta = 1 : Ntheta
    for isd = 1:Nsd
        %disp( [ 'Reading data for source ' num2str( isd ) ' of ' num2str( Nsd ) ] )
        for ird = 1:Nrcvrs_per_range
            recnum = 7 + ( itheta - 1 ) * Nsd * Nrcvrs_per_range + ( isd - 1 ) * Nrcvrs_per_range + ird - 1;
            status = fseek( fid, recnum * 4 * recl, -1 ); %Move to end of previous record
            if ( status == -1 )
                error( 'Seek to specified record failed in read_shd_bin' )
            end

            temp = fread( fid, 2 * Nrr, 'float32' );    %Read complex data
            pressure( itheta, isd, ird, : ) = temp( 1 : 2 : 2 * Nrr ) + 1i * temp( 2 : 2 : 2 * Nrr );
            %Transmission loss matrix indexed by  sd x rd x rr

        end
    end
end

fclose( fid );
