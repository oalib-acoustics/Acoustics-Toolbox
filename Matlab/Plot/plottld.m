function plottld( varargin )

% plot a single TL slice from the shade file
%
% usage:
%    plottld( filename, rrt )
%    plottld( filename, rrt, freq )
%
% where
%   filename is the shadefile (with extension)
%   rrt is the receiver range in km
%   if rrt is a vector then one plot is generated for each element
%   freq is the frequency (needed if the shdfil has multiple frequencies)
% mbp

disp( 'PlotTLd uses the first bearing and source depth in the shade file; check OK' )
itheta = 1;
isd    = 1;

narginchk( 1, 3 )

filename = varargin{1};

% optional frequency
if nargin == 2
   rrt = varargin{ 2 };
else
   rrt = NaN;
   disp( 'error: receiver range not specified' )
end

% optional frequency
if nargin == 3
   rrt  = varargin{ 2 };
   freq = varargin{ 3 };
   [ PlotTitle, ~, ~, ~, ~, Pos, pressure ] = read_shd( filename, freq );
   
else
   [ PlotTitle, ~, freqVec, ~, ~, Pos, pressure ] = read_shd( filename );
   freq = freqVec( 1 );
end


%%
% read

pressure = pressure( itheta, isd, :, : );

tlt = abs( pressure );	            % this is really the negative of TL
tlt( tlt == 0 ) = max( max( tlt ) ) / 1e10;      % replaces zero by a small number
tlt = -20.0 * log10( tlt );          % so there's no error when we take the log

% interpolate the TL field at the receiver range
% note: interp1 won't interpolate a vector with only 1 element

if ( length( Pos.r.r ) == 1 )
   tlslice = tlt;
else
   TLtemp  = squeeze( tlt )';   % need to avoid dimensional problems when TLT has 1, 2, or 3 dimensions
   tlslice = interp1( Pos.r.r, TLtemp, 1000.0 * rrt );
end

hh= plot( tlslice, Pos.r.z );
set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
set( gca, 'Xdir', 'Reverse' )
xlabel( 'TL (dB)'   )
ylabel( 'Depth (m)' )
title( { deblank( PlotTitle ); [ 'Freq = ' num2str( freq ) ' Hz    r_{rcvr} = ' num2str( rrt ) ' km' ] } )
set( hh, 'LineWidth', 2 )

% generate legend
if ( length( rrt ) > 1 )
   for irr = 1: length( rrt )
      legendstr( irr, : ) = [ 'Range = ', num2str( rrt( irr ) ), ' km' ];
   end
   
   legend( legendstr, 'Location', 'Best' )
   legend( 'boxoff' )
   drawnow
end

