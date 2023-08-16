function field( modfil, modes )

% calculates the field using modes produced by KRAKEN
% parameters read from field.flp
%
% usage: field( filename, modes )
% the filename needs to have form FileRoot.ext but the extension can be blank
% a file like 280.0_39.4 with periods in the middle will confuse it
% mbp

% Comp     = 'P';   % select component (P, H, V) (Pressure, Horizontal, Vertical)

disp( '___________________________________________________________________________' )
disp( 'Running field.m' )

clear read_modes_bin % to force rewind to beginning of mode file

[ ~, FileRoot, ext ] = fileparts( modfil );

if ( isempty( ext ) )
   ext = '.mod';   % use this as the default extension if none specified
else
   % handle cases where the modfil is FileRoot.mod.mat
   % above fileparts call only pulls off the .mat part as the extension
   if ( ext == '.mat' )
      [ ~, FileRoot, ~ ] = fileparts( FileRoot );
      ext = '.mod.mat';
   end
end

modfil = [ FileRoot ext ];

[ TitleEnv, Opt, Comp, MLimit, NProf, rProf, Pos ] = read_flp( FileRoot );

Sz   = Pos.s.z;        % source   depths
Rz   = Pos.r.z;        % receiver depths
Rr   = Pos.r.r;

NSz  = length( Sz );
NRz  = length( Rz );
NRr  = length( Rr );   % receiver ranges

freq = 0;              % assume single frequency run

% optionally read in a source beam pattern
SBP = Opt( 3 : 3 );
SrcBmPat = readpat( FileRoot, SBP );

if ( NProf == 1 )               % Range-independent case

   if nargin == 1
      Modes = read_modes( modfil, freq );
   else
      Modes = read_modes( modfil, freq, modes );
   end

   % if the title in the flp file is absent, take the title from the mode file
   if ( TitleEnv( 1 : 1 ) == '/' )
      TitleEnv = Modes.title;
   end

   freqVec = Modes.freqVec;
   Nfreq   = length( freqVec );

   MSrc    = length( Modes.k );
   M       = min( MLimit, MSrc );        % Set number of propagating modes

   pressure = zeros( Nfreq, NSz, NRz, NRr );

   % loop over frequencies
   for ifreq = 1 : Nfreq
      % following is clumsy; it has to read the mode file from the
      % beginning for each new frequency ...
      clear read_modes_bin % to force rewind to beginning of mode file

      freq = freqVec( ifreq );
      if ( ifreq > 1 )   % read new modeset
         if nargin == 1
            Modes = read_modes( modfil, freq );
         else
            Modes = read_modes( modfil, freq, modes );
         end
      end

      % calculate mode values at source depths
      % This truncates Modes.phi to exclude elastic layers
      C = interp1( Modes.z, Modes.phi( 1 : length( Modes.z ), : ), Pos.s.z );

      % loop over source depths

      for isz = 1 : NSz
         Coef = C( isz, : ).';   % modal excitation coefs for depth index isz
         if ( Modes.M == 1 )
            Coef = Coef.';   % if just one mode then Coef needs to be a row vector
         end

         % apply the source beam pattern
         if ( SBP == '*' )
            c     = 1500;   % reference sound speed, should be speed at the source depth
            omega = 2 * pi * Modes.freqVec( ifreq );
            kz2   = omega^2 / c^2 - Modes.k.^2;   % vertical wavenumber squared
            kz2( kz2 < 0 ) = 0;                   % remove negative values

            theta = atand( real( sqrt( kz2 ) ./ Modes.k ) );   % calculate the angle in degrees
            S     = interp1( SrcBmPat( :, 1 ), SrcBmPat( : , 2 ), theta );
            Coef  = Coef .* S;         % apply the shading
         end

         % calculate mode values at receiver depths
         phiR = interp1( Modes.z, Modes.phi( 1 : length( Modes.z ), : ), Pos.r.z );
         if ( Modes.M == 1 )
            phiR = phiR.';   % if just one mode then phiR needs to be a row vector
         end

         % Here's the main routine to do the modal sum
         pressure( ifreq, isz, :, : ) = evalri( Coef, phiR, Rr, Modes.k, Opt, Comp );
      end   % next source depth, isz
   end   % next frequency, ifreq
else
   if ( Opt( 2 : 2 ) == 'C' )       % Range-dependent case
      % Coupled mode theory
      clear evalcm % to force re-initialization of profile ranges
      pressure = evalcm( FileRoot, Pos, rProf, NProf, Rz, NRz, Rr, NRr, MLimit, Opt );
      freqVec( 1 ) = freq;
   else
      % Adiabatic mode theory
      error( 'Adiabatic option not implemented in Matlab version of field' )
      % pressure = evalad( FileRoot, Pos, rProf, NProf, Rz, NRz, Rr, NRr, M, Opt );
   end
end

PlotTitle = TitleEnv;
PlotType  = 'rectilin  ';
freq0     = 0.0; % Modes.freq;
atten     = 0.0;
Pos.r.r   = Rr;

save( [ FileRoot '.shd.mat' ], 'PlotTitle', 'PlotType', 'freqVec', 'freq0', 'atten', 'Pos', 'pressure' )
