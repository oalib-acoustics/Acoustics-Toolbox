function [ PlotTitle, PlotType, freqVec, freq0, atten, Pos, pressure ] = read_shd( varargin )

% Read the shade file
% calls the appropriate routine (binary, ascii, or mat file) to read in the pressure field
%
% usage: [ PlotTitle, PlotType, freqVec, freq0, atten, Pos, pressure ] = read_shd( filename );
%    Reads first source.
%        [ PlotTitle, PlotType, freqVec, freq0, atten, Pos, pressure ] = read_shd( filename, xs, ys );
%    Reads source at the specified xs, ys coordinate.
%        [ PlotTitle, PlotType, freqVec, freq0, atten, Pos, pressure ] = read_shd( filename, freq );
%    Reads source at the specified frequency.
%
% Recommended to include a file extension, if it exists.
% Otherwise it may find a different file than you intended.
%
% Output is a 5-D pressure field p( Nfreq, Ntheta, Nsd, Nrd, Nrr )
%
% If omitted, take a guess at the extension
% Matlab 'exist' command is simpler; however, it searches the whole Matlab search path.

% Determine type of file:

%error( nargchk( 1, 3, nargin, 'struct' ) );
narginchk( 1, 3 )

filename = varargin{1};

% optional frequency
if nargin == 2
    freq = varargin{ 2 };
else
    freq = NaN;
end

% optional source (x,y) coordinate
if nargin >= 3
   xs = varargin{2};
   ys = varargin{3};
else
   xs = NaN;
   ys = NaN;
end

PlotType = [];  % in case this was not set

[ ~, name, ext ] = fileparts( filename );

if ( strcmp( ext, '.mat' ) )
    [ ~, ~, ext2 ] = fileparts( name );

   switch ext2
      case '.shd'
         FileType = 'shdmat';
      case '.grn'
         FileType = 'grnmat';
   end
else
   switch filename
      case 'ASCFIL'
         FileType = 'asc';
      case 'SHDFIL'
         FileType = 'shd';
      case 'tl.grid'
         FileType = 'RAM';
      otherwise
         endchar = length( filename );
         if ( endchar >= 4 )
            FileType = lower( filename( endchar - 2 : endchar ) );
         end
   end
end

%%
switch FileType
   case { 'shd', 'grn' }   % binary format
      switch nargin
          case 1
         [ PlotTitle, PlotType, freqVec, freq0, atten, Pos, pressure ] = read_shd_bin( filename );
          case 2
         [ PlotTitle, PlotType, freqVec, freq0, atten, Pos, pressure ] = read_shd_bin( filename, freq );
          case 3
         [ PlotTitle, PlotType, freqVec, freq0, atten, Pos, pressure ] = read_shd_bin( filename, xs, ys );
      end
   case 'shdmat'   % Shade function mat file
      load( filename )
      
      % has a specific source xs, ys been given?
      if ~isnan( xs )
        xdiff = abs( Pos.s.x - xs * 1000. );
        [ ~, idxX ] = min( xdiff );
        ydiff = abs( Pos.s.y - ys * 1000. );
        [ ~, idxY ] = min( ydiff );
        
        % extract the appropriate source index
        pressureT = pressure( idxX, idxY, :, :, :, : );
        pressure = reshape( pressureT, [ Pos.Nsz, Pos.Nrz, Pos.Ntheta, Pos.Nrr ] );
      end

      % has a specific frequency been given?
      if ~isnan( freq )
         freqdiff = abs( freqVec - freq );
         [ ~, ifreq ] = min( freqdiff );
         % extract the appropriate source index
         pressure = pressure( ifreq, 1, :, : );
      end
      
   case 'asc' % ascii format
      [ PlotTitle, PlotType, freqVec, freq0, atten, Pos, pressure ] = read_shd_asc( filename );
      
   case 'grnmat'   % Green's function mat file
      load( filename )
      Pos.r.r = Pos.r.r';   % make it a column vector to match read_shd_bin

      % has a specific frequency been given?
      if ~isnan( freq )
         freqdiff = abs( freqVec - freq );
         [ ~, ifreq ] = min( freqdiff );
         % extract the appropriate source index
         pressure = pressure( ifreq, 1, :, : );
      end
      
   case 'RAM'
      [ PlotTitle, PlotType, freqVec, freq0, atten, Pos, pressure ] = read_ram_tlgrid;
   otherwise
      error( 'Unrecognized file extension' )
end

% clean up PlotTitle by taking only the part up inside the quotes
% nchars = strfind( PlotTitle, '''' );   % find quotes
% PlotTitle = [ PlotTitle( nchars( 1 ) + 1 : nchars( 2 ) - 1 ) ];
