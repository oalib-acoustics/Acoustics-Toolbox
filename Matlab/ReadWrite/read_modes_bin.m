function [ Modes ] = read_modes_bin( filename, freq, modes )

% Read the modes produced by KRAKEN
% usage:
%    [ Modes ] = read_modes_bin( filename, modes )
% filename is without the extension, which is assumed to be '.moA'
% freq is the frequency (involved for broadband runs)
%   (you can use freq=0 if there is only one frequency)
% modes is an optional vector of mode indices
%
% This version can keep reading modesets from the same modefile.
% On the first call, it opens the mode file and sets the record pointer to
% the beginning of the file.
% From then on it reads sequentially.

% derived from readKRAKEN.m    Feb 12, 1996 Aaron Thode
%
% Modes.M          number of modes
% Modes.k          wavenumbers
% Modes.z          sample depths for modes
% Modes.phi        modes
%
% Modes.Top.BC
% Modes.Top.cp
% Modes.Top.cs
% Modes.Top.rho
% Modes.Top.depth
%
% Modes.Bot.bc
% Modes.Bot.cp
% Modes.Bot.cs
% Modes.Bot.rho
% Modes.Bot.depth
%
% Modes.title      Title that was stored in the mode file
% Modes.N          Number of depth points in each medium
% Modes.Mater      Material type of each medium (acoustic or elastic)
% Modes.Nfreq      Number of frequencies
% Modes.Nmedia     Number of media
% Modes.depth      depths of interfaces
% Modes.rho        densities in each medium
% Modes.freqVec    vector of frequencies for which the modes were calculated

persistent fid iRecProfile lrecl

if ( isempty( fid ) )
   fid   = fopen( filename, 'r' );
   if ( fid == -1 )
      error( 'Mode file does not exist' )
   end
   iRecProfile = 1;   % (first time only)
   lrecl = 4 * fread( fid, 1, 'long' );   % this is converted to bytes. Fortran versions uses words instead
end

%%

rec = iRecProfile - 1;
fseek( fid, rec * lrecl + 4, -1 );

Modes.title  = fread( fid, 80, '*char' )';
Modes.Nfreq  = fread( fid,  1, 'long'  );
Modes.Nmedia = fread( fid,  1, 'long'  );
Ntot         = fread( fid,  1, 'long'  );
NMat         = fread( fid,  1, 'long'  );

if Ntot < 0, return; end

% N and Mater
rec   = iRecProfile;
fseek( fid, rec * lrecl, -1 );
for Medium = 1 : Modes.Nmedia
   Modes.N(     Medium    ) = fread( fid, 1, 'long' );
   Modes.Mater( Medium, : ) = fread( fid, 8, '*char' )';
end

%%
% depth and density

rec = iRecProfile + 1;
fseek( fid, rec * lrecl, -1 );
bulk        = fread( fid, [ 2, Modes.Nmedia ], 'float' );
Modes.depth = bulk( 1, : );
Modes.rho   = bulk( 2, : );

% frequencies
rec = iRecProfile + 2;
fseek( fid, rec * lrecl, -1 );
Modes.freqVec = fread( fid, Modes.Nfreq, 'double' );

% z
rec = iRecProfile + 3;
fseek( fid, rec * lrecl, -1 );
Modes.z = fread( fid, Ntot, 'float' );

%%
% skip through frequencies to get to the selected set

% identify the index of the frequency closest to the user-specified value
freqdiff = abs( Modes.freqVec - freq );
[ ~, freq_index ] = min( freqdiff );

% number of modes, m
iRecProfile = iRecProfile + 4;
rec = iRecProfile;

% skip through the mode file to get to the chosen frequency
for ifreq = 1 : freq_index
   fseek( fid, rec * lrecl, -1 );
   Modes.M = fread( fid, 1, 'long' );

   % advance to the next frequency
   if ( ifreq < freq_index )
      iRecProfile = iRecProfile + 3 + Modes.M + floor( 4 * ( 2 * Modes.M - 1 ) / lrecl );   % advance to next profile
      rec = iRecProfile;
   end
end

if nargin == 2
   modes = 1 : Modes.M;    % read all modes if the user didn't specify
end

% don't try to read modes that don't exist
ii    =  modes <= Modes.M;
modes = modes( ii );

%%
% Read top and bottom halfspace info

% Top
rec = iRecProfile + 1;
fseek( fid, rec * lrecl, -1 );
Modes.Top.BC    = fread( fid, 1, '*char' );

cp              = fread( fid, [ 2, 1 ], 'float' );
Modes.Top.cp    = complex( cp( 1 ), cp( 2 ) );

cs              = fread( fid, [ 2, 1 ], 'float' );
Modes.Top.cs    = complex( cs( 1 ), cs( 2 ) );

Modes.Top.rho   = fread( fid, 1, 'float' );
Modes.Top.depth = fread( fid, 1, 'float' );

% Bottom
Modes.Bot.BC    = char( fread( fid, 1, 'char' )' );

cp              = fread( fid, [ 2, 1 ], 'float' );
Modes.Bot.cp    = complex( cp( 1 ), cp( 2 ) );

cs              = fread( fid, [ 2, 1 ], 'float' );
Modes.Bot.cs    = complex( cs( 1 ), cs( 2 ) );

Modes.Bot.rho   = fread( fid, 1, 'float' );
Modes.Bot.depth = fread( fid, 1, 'float' );

%%
% Read the modes (eigenfunctions, then eigenvalues)

rec = iRecProfile;
fseek( fid, rec * lrecl, -1 );

% if there are modes, read them
if ( Modes.M == 0 )
   Modes.phi = [];
   Modes.k   = [];
else
   Modes.phi = zeros( NMat, length( modes ), 'single' );   % number of modes
   
   for ii = 1: length( modes )
      rec = iRecProfile + 1 + modes( ii );
      fseek( fid, rec * lrecl, -1 );
      phi = single( fread( fid, [ 2, NMat ], 'single' )' ); % Data is read columwise
      Modes.phi( :, ii ) = phi( :, 1 ) + 1i * phi( :, 2 );
   end
   
   % following is faster, but only valid if all the record sizes are the same
   % phitmp    = fread( fid, [ 2, Ntot * Modes.M ], 'float' ); % Data is read columwise
   % phitmp    = phitmp( 1, : ) + 1i * phitmp( 2, : );
   % Modes.phi = reshape( phitmp, Ntot, Modes.M );


   % read in the wavenumbers
   
   % Ifirst = 1;
   % cktot = [];
   
   %for I = 1 : ( 1 + ( 2 * m - 1 ) / Lrecl ),
   %   rec = 5 + m + I;
   %   fseek( fid, rec * lrecl, -1 );
   %   Ilast = min( [ m Ifirst + Lrecl / 2 - 1 ] );
   %   ck = fread( fid, [ 2, Ilast - Ifirst + 1 ], 'float' )';
   %   cktot = [ cktot; ck ];
   %   Ifirst = Ilast + 1;
   %end
   %ck = cktot( modes, 1 ) + i * cktot( modes, 2 );

   rec = iRecProfile + 2 + Modes.M;
   fseek( fid, rec * lrecl, -1 );
   
   Modes.k = zeros( 1, Modes.M );
   k       = fread( fid, [ 2, Modes.M ], 'float' );
   Modes.k = ( k( 1, : ) + 1i * k( 2, : ) ).';   % column vector
   Modes.k = Modes.k( modes );   % take the subset that the user specified
end

iRecProfile = iRecProfile + 4 + Modes.M + floor( 4 * ( 2 * Modes.M - 1 ) / lrecl );   % advance to next profile

%fclose( fid );
