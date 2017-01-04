function [ PlotTitle, PlotType, freq, atten, Pos, p ] = read_shd( filename )

% Read the shade file
% calls the appropriate routine (binary, ascii, or mat file) to read in the pressure field
% usage: [ PlotTitle, PlotType, freq, atten, Pos, p ] = read_shd( filename );
% Recommended to include a file extension, if it exists.
% Otherwise it may find a different file than you intended.
% Output is a 4-D pressure field p( Ntheta, Nsd, Nrd, Nrr )
% mbp

% if omitted, take a guess at the extension
% Matlab 'exist' command is simpler; however, it searches the whole Matlab
% search path

% Determine type of file:

PlotType = [];  % in case this was not set

switch filename
   case 'SHDFIL'
      FileType = 'shd';
   case 'SHDFIL.mat'
      FileType = 'shdmat';
   case 'ASCFIL'
      FileType = 'asc';
   case 'GRNFIL'
      FileType = 'grn';
   case 'GRNFIL.mat'
      FileType = 'grnmat';
   otherwise
      endchar = length( filename );
      if ( endchar >= 4 )
         FileType = lower( filename( endchar-2 : endchar ) );
      end
end

switch FileType
   case 'shd' % binary format
      [ PlotTitle, PlotType, freq, atten, Pos, p ] = read_shd_bin( filename );
   case 'shdmat'   % Green's function mat file
      load SHDFIL.mat
      
   case 'mat' % Matlab format
      load( filename );
   case 'asc' % ascii format
      [ PlotTitle, PlotType, freq, atten, Pos, p ] = read_shd_asc( filename );
   case 'grn' % binary format (Green's function file)
      [ PlotTitle, PlotType, freq, atten, Pos, p ] = read_shd_bin( filename );
   case 'grnmat'   % Green's function mat file
      load GRNFIL
   otherwise
      errordlg( 'Unrecognized file extension', 'PLOTSHD error' )
      return
end

% clean up PlotTitle by taking only the part up inside the quotes
% nchars = strfind( PlotTitle, '''' );   % find quotes
% PlotTitle = [ PlotTitle( nchars( 1 ) + 1 : nchars( 2 ) - 1 ) ];