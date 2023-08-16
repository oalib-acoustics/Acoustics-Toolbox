function fieldsco( filename, PlotTitle, freq, atten, Pos, Gtemp, Opt, Rminkm, Rmaxkm, Nrr )

% Calculates the field using the Green's function produced by SCOOTER
%
% This version uses the trapezoidal rule directly to do a DFT, rather than an FFT
%
% You could save storage by doing the transform for one row at a time out
% of the Green's function file, i.e. a single source/receiver depth.
% This version does the entire matrix at once. If that fits comfortably in
% memory then it is potentially faster in terms of using multiple cores
% efficiently
%
% usage:
%    fieldsco( filename )
%  or
%    fieldsco( filename, PlotTitle, freq, atten, Pos, Gtemp, option, Rminkm, Rmaxkm, Nrr )
%
% You must include the file extension, if it exists
%
% Optionally reads control info from a file fields.flp
% mbp, Dec. 2012 based on fieldsco

% phase speed limits for optional tapering
% this is for user play (at your own risk)

global cmin cmax

if ( isempty( cmin ) )
   cmin = 1e-10;
   cmax = 1e30;
end

% if user is just supplying a file name then read all the variables

if ( strcmp( filename( end - 3 : end ), '.mat' ) )
   fileroot = filename( 1 : end - 8 );
else
   fileroot = filename( 1 : end - 4 );
end

% disp( [ 'fieldsco ' filename ] )

if ( nargin == 1 )
   % read fields.flp data to select range limits
   fid = fopen( [ fileroot '.flp' ], 'r' );
   if ( fid < 3 )
      error( 'flp file missing' )
   end
   
   Opt = fgetl( fid );   % x or r coordinate; positive/negative/both sides of spectrum
   % Extract letters between the quotes
   nchars = strfind( Opt, '''' );   % find quotes
   Opt    = Opt( nchars( 1 ) + 1 : nchars( 2 ) - 1 );
   if ( length( Opt ) <= 3 )
      Opt( 4 : 4 ) = 'O';   % default beampattern is omni
   end

   disp( 'Option line:' )
   disp( Opt )

   % read the range vector
   %[ rr_km, Nrr ] = readvector( fid );
   Rr = readr( fid );
   Nrr = length( Rr );
   Rr = 1000 * Rr;   % convert km to m
   Rr( abs( Rr ) < realmin ) = 1;   % replace zero range with 1 m to avoid singularity
   % fprintf( '\n\n--- FieldSCO --- \n\nRminkm = %d, Rmaxkm = %d, Nrr = %i \n', Rminkm, Rmaxkm, Nrr )
   
   fclose( fid );
   
   % read in the Green's function; just getting basic info at this stage
   [ ~, ~, freqVec, freq0, ~, Pos, ~ ] = read_shd( filename );
end

% optionally read in a source beam pattern
SBP = Opt( 4 : 4 );
SrcBmPat = readpat( fileroot, SBP );

Nsz   = length( Pos.s.z );    % # of source   depths
Nrz   = length( Pos.r.z );    % # of receiver depths
Nfreq = length( freqVec );    % # of frequencies

pressure = zeros( Nfreq, Nsz, Nrz, Nrr );

for ifreq = 1 : Nfreq
   freq   = freqVec( ifreq );
   omega  = 2.0 * pi * freq;
   kleft  = omega / cmax; % left  limit for tapering
   kright = omega / cmin; % right limit for tapering
   
   % read in the Green's function for this frequency
   [ PlotTitle, ~, freqVec, freq0, atten, Pos, Gtemp ] = read_shd( filename, freqVec( ifreq ) );

   cVec = Pos.r.r;   % the phase speeds for the wavenumbers are stored in Pos.r.r
   k    = 2 * pi * freq  ./ cVec;

   % for a SPARC run the wavenumber vector is independent of frequency
   if ( contains( PlotTitle( 1 : 6 ), 'SPARC' ) )
      k    = 2 * pi * freq0 ./ cVec;
   end
   
   Nk   = length( k );              % # of wavenumbers
   
   deltak = ( k( end ) - k( 1 ) ) / ( length( k ) - 1 );

   if ( Rr( end ) * deltak > 10 )
      error( 'The wavenumber sampling is too coarse to accurately calculate the field at the largest range' )
   end
   
   if ( contains( PlotTitle( 1 : 7 ), 'SCOOTER' ) )
      atten  = deltak;   % this is calculated here because SCOOTER doesn't write it for each frequency
   end
   
   ck = k + 1i * atten;
   x  = ck.' * abs( Rr' );   % this is a matrix of size nk * nrr
   
   switch Opt( 1 : 1 )
      case 'X'   % line source
         X  = exp( -1i * x );   % e^( -i k r ) matrix
         X2 = exp( +1i * x );   % e^( +i k r ) matrix
         
         factor1 = 1;
         factor2 = deltak ./ sqrt( 2 * pi );
      case 'R'   % point source
         X  = exp( -1i * ( x - pi / 4 ) );   % e^( -i k r ) matrix
         X2 = exp( +1i * ( x - pi / 4 ) );   % e^( +i k r ) matrix
         
         factor1 = sqrt( ck );
         factor2 = deltak ./ sqrt( 2 * pi * abs( Rr' ) );
      case 'S'   % point source with cylindrical spreading removed
         X  = exp( -1i * ( x - pi / 4 ) );   % e^( -i k r ) matrix
         X2 = exp( +1i * ( x - pi / 4 ) );   % e^( +i k r ) matrix
         
         factor1 = sqrt( ck );
         factor2 = deltak ./ sqrt( 2 * pi );
      case 'H'   % exact Bessel transform
         X  = besselh( 0, 1, x );   % Hankel function
         X2 = besselh( 0, 2, x );   % Hankel function
         
         factor1 = ck;
         factor2 = deltak / 2;
      otherwise
         disp( 'fieldsco.m -- Unknown first letter option in second line of .flp file' )
   end
   
   for isz = 1 : Nsz
      G = squeeze( Gtemp( isz, 1, :, : ) );

      if size( G, 2 ) == 1 % if G is a vector, expand it into a matrix with one row
         G = reshape( G, 1, length( G ) );
      end
      
      % apply the source beam pattern
      if ( SBP == '*' )
         c = 1500;   % reference sound speed, should be speed at the source depth
         kz2 = omega^2 / c^2 - k.^2;   % vertical wavenumber squared
         kz2( kz2 < 0 ) = 0;           % remove negative values
         
         theta = atand( sqrt( kz2 ) ./ k );   % calculate the angle in degrees
         S = interp1( SrcBmPat( :, 1 ), SrcBmPat( : , 2 ), theta );
         G = scalecol( G, S );         % apply the shading
      end
      
      G = taper( G, k, Nk, kleft, kright );
      G( abs( G ) < 1e-40 ) = 0;   % avoid underflows--- they slow the following down a huge amount
      
      switch Opt( 2: 2 )
         case 'P'
            G = scalecol( G, factor1 );
            Y = -G * X;   % here's where the DFT is done
         case 'N'   % check this !!!
            G = scalecol( G, factor1 );
            Y = -G * X2;
         case 'B'   % check this !!!
            G = scalecol( G, factor1 );
            Y = -G * ( X + X2 );
            %Y = -G * cos(( x - pi / 4 ) );
         otherwise
            disp( 'fieldsco.m -- Unknown second letter option in second line of .flp file' )
      end
      
      pressure( ifreq, isz, :, : ) = scalecol( Y, factor2 );  % cylindrical spreading, etc.
     
      %figure; imagesc( rr, Pos.r.z, 20 * log10( abs( squeeze( pressure( ifreq, 1, :, : ) ) ) ) )
%       figure; imagesc( rr, Pos.r.z, abs( squeeze( pressure( ifreq, 1, :, : ) ) ) )
%       title( num2str( freqVec( ifreq ) ) )
      fprintf( '      Transform completed for source depth %f, frequency (or time) %f \n', Pos.s.z( isz ), freq );
   end   % next source depth
end   % next frequency

Pos.r.r  = Rr;
atten    = 0;
PlotType = 'rectilin  ';

save( [ fileroot '.shd.mat' ], 'PlotTitle', 'PlotType', 'freqVec', 'freq0', 'atten', 'Pos', 'pressure' )


%%

function G = taper( G, k, Nk, kleft, kright )

% windowing to smooth out any discontinuities at the end of the spectrum

if ( kleft > k( 1 )  )
   Nwinleft = 2 * round( ( kleft - k( 1 ) ) / ( k(end) - k( 1 ) ) * Nk ) + 1;  % odd number covering 10% of the spectrum
   winleft  = hanning( Nwinleft )';
   Nwinleft = ( Nwinleft - 1 ) / 2;
   winleft  = winleft( 1 : Nwinleft ); % taking just the left half of the symmetric Hanning window
else
   Nwinleft = 0;
   winleft  = [];
end

if ( kright < k( end )  )
   Nwinright = 2 * round( ( k(end) - kright ) / ( k(end) - k( 1 ) ) * Nk ) + 1; % odd number covering 10% of the spectrum
   winright  = hanning( Nwinright )';
   Nwinright = ( Nwinright - 1 ) / 2;
   winright  = winright( end - Nwinright + 1 : end ); % taking just the right half of the symmetric Hanning window
else
   Nwinright = 0;
   winright  = [];
end

if ( Nk < Nwinleft + Nwinright )
   error( [ mfilename, 'phase speed limits for windowing are invalid' ] );
end

window  = [ winleft ones( 1, Nk - Nwinleft - Nwinright ) winright ];

G = scalecol( G, window );
