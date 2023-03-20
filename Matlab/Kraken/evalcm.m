function P = evalcm( FileRoot, Pos, RProf, NProf, rd, Nrd, r, Nr, M, Opt )

% conversion from Fortran routine evalcm.f90
% mbp 4/2009


% Computes pressure field using coupled mode theory
% Normalized to pressure of point source at 1 meter

% Opt = X     Cartesian   (X, z) coordinates
% Opt = R     Cylindrical (R, z) coordinates

% Note number of propagating modes is reset after first segment.
% Thus M restricts the number of modes in the source field but
% thereafter energy is allowed to couple into higher-order modes.

% Should half-space contribution include imaginary part?
% Then we need to handle the growing part of leaky waves.

% Note: Later work with a PE suggests that it may be important to do a
% spline fit for the projection of the pressure at each new segment

persistent first
clear read_modes_bin % to force re-initialization of record counter in mode file

P = zeros( 1, 1, Nrd, Nr );   % pre-allocated pressure field

% Compute ranges (in meters) where new profiles are used

if ( isempty( first ) )
   for iProf = NProf : -1 : 2
      RProf( iProf ) = 500.0 * ( RProf( iProf ) + RProf( iProf - 1 ) );
   end
   
   RProf(  NProf + 1 ) = realmax;
   first = 0;
end

% Evaluate mode excitation coefficients, A(mode)
iProf    = 1;
filename = [ FileRoot '.mod' ];
Modes    = read_modes( filename, 0.0 );   % freq=0.0 because broadband is not allowed here

MSrc     = length( Modes.k );

% weight modes by mode excitation
'mode sampling should be done using interp1 !!!'

zs   = Pos.s.z;
isd  = find( Modes.z >= zs );    % index of source depth
isd  = isd( 1 );
phiS = Modes.phi( isd, : ).';    % column vector with modal weights

%M = min( MLimit, MSrc );        % Set number of propagating modes
k     = Modes.k;
phi   = interp1( Modes.z, Modes.phi, rd' );   % subsample phi onto the receiver grid
M1    = length( k );
M     = min( M, M1 );
k     = k( 1 : M );   % reduce size of vector based on Mlimit

% Initialize the amplitudes A
if ( Opt(1:1) == 'X' )   % Cartesian coordinates
   A =      sqrt( 2.0 * pi ) * exp( 1i * pi / 4.0 ) * phiS( 1:M ) ./       k( 1:M );
else                     % Cylindrical coordinates
   A = 1i * sqrt( 2.0 * pi ) * exp( 1i * pi / 4.0 ) * phiS( 1:M ) ./ sqrt( k( 1:M ) );
end

%A = zeros( length( A ), 1 );
%A( 500 ) = 1;

% March forward in range
for ir = 1 : Nr
   if ( r( ir ) > RProf( iProf + 1 ) )  % Crossing into new range segment?
      
      iProf = iProf + 1;
      
      % Advance to interface
      if ( ir == 1 )
         A = A .* exp( -1i * k * RProf( iProf ) );                   % advance the phase of the coefficients
      else
         A = A .* exp( -1i * k * ( RProf( iProf ) - r( ir - 1 ) ) ); % advance the phase of the coefficients
      end
      
      % Here's where we cross over
      if ( iProf  <= NProf )
         [ Modes, A, k, phi, M ] = NewProfile( FileRoot, iProf, Modes, M, rd, A );
         disp( 'New profile read' )
         disp( [ r(ir ), iProf ] )
         %disp( ' #modes = ' )
         %disp( M )
      end
      
      % Are there other segments to cross?
      while ( r( ir ) > RProf( iProf + 1 ) )
         iProf = iProf + 1;
         A = A .* exp( -1i * k * ( RProf( iProf ) - RProf( iProf - 1 ) ) ); % advance the phase of the coefs.
         if ( iProf <= NProf )
            [ Modes, A, k, phi, M ] = NewProfile( FileRoot, iProf, Modes, M, rd, A );
         end
      end
      
      % Advance the remaining distance
      A = A .* exp( -1i * k * ( r( ir ) - RProf( iProf ) ) );     % advance the phase of the coefficients
      
   else

      if ( ir == 1 )
         A = A .* exp( -1i * k * r( ir ) );                       % advance the phase of the coefficients
      else
         A = A .* exp( -1i * k * ( r( ir ) - r( ir - 1 ) ) );     % advance the phase of the coefficients
      end
      
   end

   % Add up modal contributions
   if ( M == 0 )    % modes cut off?
      P( 1, 1, :, ir ) = 0;
   else
      if ( Opt(1:1) == 'R' && r( ir ) ~= 0.0 )
         P( 1, 1, :, ir ) =  phi( :, 1 : M ) * A( 1 : M ) / realsqrt( r( ir ) );
      else
         P( 1, 1, :, ir ) =  phi( :, 1 : M ) * A( 1 : M );
      end
   end

end    % next range step

%**********************************************************************C
function [ Modes, A, kR, phiR, MR ] = NewProfile( FileRoot, iProf, ModesL, M, rd, A )

% For a given profil number:
%     read in modes for current segment
%     project the pressure field onto the new modes
%     extract values of the modes at rcvr depths

% Compute pressure along the left of the interface and read new mode set

if ( M == 0 ) % quick return if no propagating modes
   Modes = [];
   A     = [];
   kR    = [];
   phiR  = [];
   MR    = 0;
   return
end

[ P, z, kR, NR, NTot, ML, MR, gamTL, gamBL, depthTL, depthBL, phiTL, phiBL, Modes ] = ...
   pleft( FileRoot, iProf, ModesL, A, M );

% Calculate new amplitudes by projecting pressure onto the new modes
if ( MR == 0 ) % quick return if no propagating modes
   phiR = [];
   kR   = [];
   return
end

A = zeros( MR, 1 );   % shrink size of A as necessary

% calculate modal values for top and bottom halfspaces
if ( Modes.Top.BC == 'A' )   % Top halfspace
   phiTR = Modes.Top.phi;
   gamTR = Modes.Top.gamma;
end

if ( Modes.Bot.BC == 'A' )   % Bottom halfspace
   phiBR = Modes.Bot.phi;
   gamBR = Modes.Bot.gamma;
end

for mode = 1 : MR
   % tabulate the new mode on the grid from the previous segment
   phiTmp = Modes.phi( :, mode );
   
   if ( Modes.Bot.BC == 'A' )
      ii = find ( z > Modes.Bot.depth );
      phiTmp( ii ) = phiBR( mode ) * exp( -gamBR( mode ) * ( z( ii ) - Modes.Bot.depth ) );
   end
   
   if ( Modes.Top.BC == 'A' )
      ii = find ( z < Modes.Top.depth );
      phiTmp( ii ) = phiTR( mode ) * exp( -gamTR( mode ) * ( Modes.Top.depth - z( ii ) ) );
   end
   
   % Compute new amplitudes:
   %    A = Integral[ P( z ) * phi( z ) dz ]
   %    Integral is done using trapezoidal rule
   
   sum1 = P.' * phiTmp;
   
   if ( Modes.Top.BC == 'A' )  % contribution from upper halfspace
      sum1   = sum1 + tail( z( 1    ), phiTL, gamTL, depthTL, ...
         phiTR( mode ) / Modes.Top.rho, gamTR( mode ), Modes.Top.depth );
   end
   
   if ( Modes.Bot.BC == 'A' )  % contribution from lower halfspace
      sum1   = sum1 + tail( z( NTot ), phiBL, gamBL, depthBL, ...
         phiBR( mode ) / Modes.Bot.rho, gamBR( mode ), Modes.Bot.depth );
   end
   A( mode ) = sum1;
end    % next mode

% Subtabulate modes at receiver depths

phiR( :, : ) = interp1q( Modes.z, Modes.phi( :, 1 : MR ), rd' );   % subsample phi onto the receiver grid
if ( MR == 1 )   % interp1 converted phiR to a row vector if MR == 1
   phiR = phiR.';
end

if ( Modes.Top.BC == 'A' )
   ii = find ( rd < Modes.Top.depth ); % Rcvr in upper halfspace
   phiR( ii, : ) = phiTR .* exp( -gamTR * ( Modes.Top.depth   - rd( ii  ) ) ).';
end

if ( Modes.Bot.BC == 'A' )
   ii = find ( rd > Modes.Bot.depth ); % Rcvr in lower halfspace
   phiR( ii, : ) = ( diag( phiBR ) * exp( -gamBR * ( rd( ii )  - Modes.Bot.depth   ) ) ).';
end


fprintf( 'depth-averaged power: %10.2e \n', norm( A( 1:MR ) )^2 )

%**********************************************************************C

function [ P, z, kR, NR, NTot, ML, M, gamTL, gamBL, depthTL, depthBL, phiTL, phiBL, Modes ] = ...
   pleft( FileRoot, iProf, Modes, A, ML )

% Computes the pressure field along the interface using the depth sampling of the
% next segment.
% Also returns information needed for the tails in the halfspaces

% Get modal info at end of last segment (stored in structure Modes)

zL      = Modes.z;
depthBL = Modes.Bot.depth;
depthTL = Modes.Top.depth;
ML      = min( length( Modes.k ), ML );   % we only know M amplitudes
A       = A( 1 : ML );
Modes.k = Modes.k( 1 : ML );
NL      = length( Modes.z );

% Compute pressure at the interface

% Halfspace information
phiTL = A .* Modes.phi( 1,  1 : ML ).';
phiBL = A .* Modes.phi( NL, 1 : ML ).';

gamTL = zeros( ML, 1 );
gamBL = zeros( ML, 1 );

if ( Modes.Top.BC == 'A' )   % Top halfspace
   gamTL  = Modes.Top.gamma;
end

if ( Modes.Bot.BC == 'A' )   % Bottom halfspace
   gamBL  = Modes.Bot.gamma;
end

PL =  Modes.phi( :, 1 : ML ) * A;

% Read modal data in new segment
filename = [ FileRoot '.mod' ];
Modes = read_modes( filename, 0.0 );   % freq=0.0 because broadband is not allowed here

z  = Modes.z;

NR = length( z );
M  = Modes.M;

if ( M == 0 )   % jump out if modefile was empty
   kR   = [];
   NR   = [];
   NTot = [];
   P    = zeros( size( z ) );
   return
end

% Upslope? Extend the z vector with data from zL
NTot = NR;
ii   = find( zL > Modes.z( NTot ) );
newz = length( ii );

z( NTot + 1 : NTot + newz ) = zL( ii );
NTot = NTot + newz;

kR = Modes.k( 1 : M );
if ( Modes.z( 1 ) ~= Modes.Top.depth || Modes.z( end ) ~= Modes.Bot.depth )
   disp( 'Fatal Error: modes must be tabulated throughout the ocean and sediment to compute the coupling coefs.' )
   stop
end

% Tabulate medium density on new grid

rhomed = zeros( size( z ) );

for med = 1 : Modes.Nmedia - 1
   ii = find( z >= Modes.depth( med ) & z <= Modes.depth( med + 1 ) );
   rhomed( ii ) = Modes.rho( med );
end

% special case for last medium
ii = find( z >= Modes.depth( Modes.Nmedia ) );
rhomed( ii ) = Modes.rho( Modes.Nmedia );

% special case for upper halfspace
ii = find( z < Modes.Top.depth );
rhomed( ii ) = Modes.Top.rho;

% special case for lower halfspace
ii = find( z > Modes.Bot.depth );
rhomed( ii ) = Modes.Bot.rho;

% interpolate pressure onto the new grid
P = interp1q( zL, PL, z );

% special case for upper halfspace
if ( Modes.Top.BC == 'A' )
   ii = find( z < depthTL );
   if ( ~isempty( ii ) )
      P( ii ) = phiTL( 1 : ML ).' * exp( -gamTL( 1 : ML ) * ( depthTL - z( ii )' ) );
   end
end

% special case for lower halfspace
if ( Modes.Bot.BC == 'A' )
   ii = find( z > depthBL );
   if ( ~isempty( ii ) )
      P( ii ) = phiBL( 1 : ML ).' * exp( -gamBL( 1 : ML ) * ( z( ii )' - depthBL ) );
   end
end

% compute mesh width, h (it's actually h/rho)

% check this is correct for an irregular mesh

h = zeros( size( z ) );

h( 1           ) = 0.5 * ( z( 2       ) - z( 1           ) )  / rhomed( 1           ); % forward  difference
h( 2 : end - 1 ) = 0.5 * ( z( 3 : end ) - z( 1 : end - 2 ) ) ./ rhomed( 2 : end - 1 ); % centered difference
h( NTot        ) = 0.5 * ( z( NTot    ) - z( NTot - 1    ) )  / rhomed( NTot        ); % backward difference

% Point just above or below the interface
for med = 2 : Modes.Nmedia
   % find points where the centered difference straddles an interface
   % typically 3 points will give two straddlers
   % hence the ./ rather than / below
   % see my notes
   DBelow  = Modes.depth( med );
   ii = find( z( 1 : end - 2 ) < DBelow & z( 3 : end ) > DBelow ) + 1;
   
   h( ii ) = 0.5 * ( z( ii + 1 ) ./ rhomed( ii + 1 ) - z( ii - 1 ) ./ rhomed( ii - 1 ) ...
                        - DBelow ./ rhomed( ii + 1 ) +      DBelow ./ rhomed( ii - 1 ) );
end

P = h .* P;   % weight P by that h for later use in Trapezoidal rule


%**********************************************************************C
function [ tail ] = tail( D, phiL, gamL, DL, phiR, gamR, DR )

% Computes the contributions from the tails in the halfspaces
% gamL and phiL are vectors

FR =  phiR * exp( -gamR * ( D - DR ) );

if ( D == DL )
   tail = FR * sum( phiL                              ./ ( gamL + gamR ) );
else
   tail = FR * sum( phiL .* exp( -gamL * ( D - DL ) ) ./ ( gamL + gamR ) );
end

function rootz = PekerisRoot( z )

% Calculates the Pekeris branch of the square root
ii = find( real( z ) >= 0.0 );
rootz( ii ) = sqrt( z( ii ) );

ii = find( real( z ) < 0.0 );
rootz( ii ) = 1i * sqrt( -z( ii ) );

rootz = rootz( : );   % make output a column vector

