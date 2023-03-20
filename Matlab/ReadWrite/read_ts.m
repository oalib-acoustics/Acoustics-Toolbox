function [ PlotTitle, Pos, tout, RTS ] = read_ts( filename )

% Read the time-series file
% calls the appropriate routine (binary, ascii, or mat file) to read in the pressure field
%
% usage: [ PlotTitle, PlotType, freqVec, atten, Pos, pressure ] = read_shd( filename );
%    Reads first source.
%

% Recommended to include a file extension, if it exists.
% Otherwise it may find a different file than you intended.
%
% If omitted, take a guess at the extension
% Matlab 'exist' command is simpler; however, it searches the whole Matlab search path.

% Determine type of file:

[ ~, ~, ext ] = fileparts( filename );

if ( strcmp( ext, '.mat' ) )
   load( filename, 'PlotTitle', 'Pos', 'tout', 'RTS' )
   RTS = RTS';
else
   % open the file
   fid = fopen( filename, 'r' );
   if ( fid == -1 )
      error( 'No timeseries file with that name exists' );
   end
   
   % read
   
   PlotTitle = fgetl( fid );
   nrd      = fscanf( fid, '%f', 1 );
   rd       = fscanf( fid, '%f', nrd );
   temp     = fscanf( fid, '%f', [ nrd + 1, inf ] );
   fclose( fid );
   
   % extract rts
   Pos.r.z = rd;
   tout = temp( 1, : )';
   RTS  = temp( 2 : nrd + 1, : )';
end
