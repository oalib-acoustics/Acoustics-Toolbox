function P = evalcmLoadAll( fileroot, Pos, Dwedge, rd, Nrd, r, Nr, M, Opt )

% conversion from Fortran routine evalcm.f90
% mbp 4/2009

% Computes pressure field using coupled mode theory
% Normalized to pressure of point source at 1 meter

% Opt = X     Cartesian   (X, z) coordinates
% Opt = R     Cylindrical (R, z) coordinates

% Number of propagating modes is reset after first segment.
% Thus M restricts the number of modes in the source field but
% thereafter energy is allowed to couple into higher-order modes.

% This version reads all the modes in at the beginning, then reuses as
% needed based on the water depth

% Dwedge defines the depths associated with each mode file
% The ith modefile must correspond to Dwedge( i )

% Note: Later work with a PE suggests that it may be important to do a
% spline fit for the projection of the pressure at each new segment

global xBot

P = zeros( Nrd, Nr );   % pre-allocate pressure field

% ********* load all the mode information
MaxM = 2000
phiR = zeros( Nrd, MaxM, length( Dwedge ), 'single' );

for ii = 1 : length( Dwedge )
    filename = sprintf( '%s_rd.mod', fileroot );
    fprintf( 'Reading mode file # %4i %s \r', ii, filename );
    
    Modes( ii ) = read_modes( filename, 0.0 );   % dummy frequency=0.0
    
    % Subtabulate modes at receiver depths
    % Values outside the boundaries are set to zero (except for halfspaces)
    phiR( :, 1:size( Modes( ii ).phi, 2 ), ii ) = ...
        interp1( Modes( ii ).z, Modes( ii ).phi, rd', 'linear', 0.0 );
    
    %    if ( MR == 1 )   % interp1 converted phiR to a row vector if MR == 1
    %       phiR = phiR.';
    %    end
    
     iTop = find ( rd < Modes( ii ).Top.depth ); % Rcvr in upper halfspace
     phiR( iTop, 1 : Modes( ii ).M, ii ) = exp( -( Modes( ii ).Top.depth   - rd( iTop  ) ).' ...
                   * Modes( ii ).Top.gamma.' ) * diag( Modes( ii ).Top.phi );
    
     iBot = find ( rd > Modes( ii ).Bot.depth ); % Rcvr in lower halfspace
     phiR( iBot, 1 : Modes( ii ).M, ii ) = exp( -( rd( iBot )  - Modes( ii ).Bot.depth   ).' ...
                   * Modes( ii ).Bot.gamma.' ) * diag( Modes( ii ).Bot.phi );
end

% get modes at first range
DOld = interp1q( xBot( 1, : )', xBot( 2, : )', r( 1 ) );

% Find first profile depth <= to the depth at the source
% assumes upslope ording of wedge modes
% assumes the wedge modes cover all possible depths
iDNew = find ( Dwedge <= DOld );

iDOld = iDNew( 1 );
k     = Modes( iDOld ).k;
M1    = length( k );
M     = min( M, M1 );
k     = k( 1 : M );   % reduce size of vector based on Mlimit

% evaluate modes at source depth
zs   = Pos.s.z;
phiS = interp1q( Modes( ii ).z, Modes( iDOld ).phi, zs ).';   % subsample phi onto the receiver grid

% Initialize the amplitudes A
if ( Opt(1:1) == 'X' )   % Cartesian coordinates
    A =      sqrt( 2.0 * pi ) * exp( 1i * pi / 4.0 ) * phiS( 1 : M ) ./       k( 1 : M );
else                     % Cylindrical coordinates
    A = 1i * sqrt( 2.0 * pi ) * exp( 1i * pi / 4.0 ) * phiS( 1 : M ) ./ sqrt( k( 1 : M ) );
end

% March forward in range
for ir = 1 : Nr
    % calculate the bottom depth at this range
    DNew  = interp1q( xBot( 1, : )', xBot( 2, : )', r( ir ) );
    
    % Indices of sections between the depth limits
    iDNew = find ( Dwedge >= min( DOld, DNew ) & Dwedge <= max( DOld, DNew ) );
    
    % Order results, upslope vs. downslope)
    if ( DNew > DOld )   % going downslope?
        iDNew = fliplr( iDNew );
    end
     
    if ( ir == 1 )   % first receiver range is a special case
        phi = phiR( :, 1 : M, iDOld );          % get modes at receiver depths
        A   = A .* exp( -1i * k * r( ir ) );  % advance the phase of the coefficients
        
    else
        if ( ~isempty( iDNew ) && iDNew( end ) ~= iDOld )   % Crossing into new profile?
            
            % step through each interface
            rOld  = r( ir - 1 );
            DOld0 = Modes( iDOld ).Bot.depth;

            for istep = 1 : length( iDNew )
                % calculate range of the interface
                alpha = ( Modes( iDNew( istep ) ).Bot.depth - DOld0 ) / ...
                        ( DNew                              - DOld0 ); % proportional distance
                rNew  = ( 1 - alpha ) * r( ir - 1 ) + alpha * r( ir ); % range of the interface
                
                % disp( [ rNew DOld0 DNew Modes( iDNew( istep ) ).Bot.depth alpha ] )
                
                A     = A .* exp( -1i * k * ( rNew - rOld ) );    % Advance to the interface
                [ A, k, M ] = NewProfile( Modes( iDOld ), Modes( iDNew( istep ) ), M, A );
                
                rOld  = rNew;
                iDOld = iDNew( istep );
            end
            
            A     = A .* exp( -1i * k * ( r( ir ) - rOld ) );   % Advance the remaining distance
            DOld  = DNew;
            iDOld = iDNew( end );
            phi   = phiR( :, 1:M, iDOld );                      % get new modes at receiver depths
        else
            A = A .* exp( -1i * k * ( r( ir ) - r( ir - 1 ) ) );  % advance the phase of the coefficients
        end
        
    end
    
    % Add up modal contributions
    if ( M == 0 )    % modes cut off?
        P( 1, 1, :, ir ) = 0;
    else
        if ( Opt(1:1) == 'R' && r( ir ) ~= 0.0 )
            P( :, ir ) =  phi * A / realsqrt( r( ir ) );
        else
            P( :, ir ) =  phi * A( 1 : M );
        end
    end
    
end    % next range step

P = reshape( P, 1, 1, Nrd, Nr );


%**********************************************************************
function [ A, kR, MR ] = NewProfile( ModesL, ModesR, M, A )

% For a given profil number:
%     read in modes for current segment
%     project the pressure field onto the new modes
%     extract values of the modes at rcvr depths

% Compute pressure along the left of the interface and read new mode set

if ( M == 0 ) % quick return if no propagating modes
    A     = [];
    kR    = [];
    MR    = 0;
    return
end

[ P, z, kR, NTot, MR, gamTL, gamBL, depthTL, depthBL, phiTL, phiBL ] = ...
    pleft( ModesL, ModesR, A, M );

% Calculate new amplitudes by projecting pressure onto the new modes
if ( MR == 0 ) % quick return if no propagating modes
    return
end

A = zeros( MR, 1 );   % shrink size of A as necessary

% calculate modal values for top and bottom halfspaces

iTop = find ( z < ModesR.Top.depth );   % indices of depths in Top halfspace
iBot = find ( z > ModesR.Bot.depth );   % indices of depths in Bot halfspace
phiTR = ModesR.Top.phi;
gamTR = ModesR.Top.gamma;
phiBR = ModesR.Bot.phi;
gamBR = ModesR.Bot.gamma;

sum.Top.Tails = zeros( MR, 1 );
if ( ModesR.Top.bc == 'A' )   % Top halfspace
    sum.Top.Tails   = tail( z( 1    ), phiTL, gamTL, depthTL, ...
        phiTR( mode ) / ModesR.Top.rho, gamTR( mode ), ModesR.Top.depth );
end

sum.Bot.Tails = zeros( MR, 1 );
if ( ModesR.Bot.bc == 'A' )   % Bottom halfspace
    sum.Bot.Tails   = tail( z( NTot ), phiBL, gamBL, depthBL, ...
        phiBR / ModesR.Bot.rho, gamBR, ModesR.Bot.depth );
end

% Main loop: project pressure field on each mode in the new segment

for imode = 1 : MR
    % tabulate the new mode on the grid from the previous segment
    phi = ModesR.phi( :, imode );
    
    % following should be pulled out of the loop for efficiency?
    phi( iTop ) = 0.;
    phi( iBot ) = 0.;

    if ( ModesR.Top.bc == 'A' )
       phi( iTop ) = phiTR( imode ) * exp( -gamTR( imode ) * ( ModesR.Top.depth - z( iTop ) ) );
    end
    
    if ( ModesR.Bot.bc == 'A' )
       phi( iBot ) = phiBR( imode ) * exp( -gamBR( imode ) * ( z( iBot ) - ModesR.Bot.depth ) );
    end
    
    % Compute new amplitudes:
    %    A = Integral[ P( z ) * phi( z ) dz ]
    %    Integral is done using trapezoidal rule
    
    A( imode ) = P.' * phi + sum.Top.Tails( imode ) + sum.Bot.Tails( imode );
    
end    % next mode

fprintf( 'depth-averaged power: %10.2e \n', norm( A( 1:MR ) )^2 )

%**********************************************************************

function [ P, z, kR, NTot, M, gamTL, gamBL, depthTL, depthBL, phiTL, phiBL ] = ...
    pleft( ModesL, ModesR, A, ML )

% Computes the pressure field along the interface using the depth sampling of the
% next segment.
% Also returns information needed for the tails in the halfspaces

% Get modal info at end of last segment (stored in structure Modes)

zL      = ModesL.z;
depthBL = ModesL.Bot.depth;
depthTL = ModesL.Top.depth;
ML      = min( length( ModesL.k ), ML );   % we only know M amplitudes
A       = A( 1 : ML );
Modes.k = ModesL.k( 1 : ML );
NL      = length( ModesL.z );

% Halfspace information
phiTL = A.' .* ModesL.phi( 1,  1 : ML );
phiBL = A.' .* ModesL.phi( NL, 1 : ML );

gamTL = zeros( 1, ML );
gamBL = zeros( 1, ML );

if ( ModesL.Top.bc == 'A' )   % Top halfspace
    gamTL  = ModesL.Top.gamma;
end

if ( ModesL.Bot.bc == 'A' )   % Bottom halfspace
    gamBL  = ModesL.Bot.gamma;
end

% phiTemp = ModesL.phi( :, 1 : ML );
PL = ModesL.phi( :, 1 : ML ) * A;   % Compute pressure at the interface

z  = ModesR.z;

NR = length( z );
M  = ModesR.M;

if ( M == 0 )   % jump out if modefile was empty
    kR   = [];
    NTot = [];
    P    = zeros( size( z ) );
    return
end

% Upslope? Extend the z vector with data from zL
NTot = NR;

ii   = find( zL > ModesR.z( NTot ) );
newz = length( ii );

z( NTot + 1 : NTot + newz ) = zL( ii );
NTot = NTot + newz;

kR = ModesR.k( 1 : M );
if ( ModesR.z( 1 ) ~= ModesR.Top.depth || ModesR.z( end ) ~= ModesR.Bot.depth )
    disp( 'Fatal Error: modes must be tabulated throughout the ocean/sediment to compute the coupling' )
    stop
end

% Tabulate medium density on new grid

rhomed = zeros( size( z ) );

for med = 1 : ModesR.Nmedia - 1
    ii =  z >= ModesR.depth( med ) & z <= ModesR.depth( med + 1 ) ;
    rhomed( ii ) = ModesR.rho( med );
end

% special case for last medium
ii =  z >= ModesR.depth( ModesR.Nmedia ) ;
rhomed( ii ) = ModesR.rho( ModesR.Nmedia );

% special case for upper halfspace
ii =  z < ModesR.Top.depth ;
rhomed( ii ) = ModesR.Top.rho;

% special case for lower halfspace
ii =  z > ModesR.Bot.depth ;
rhomed( ii ) = ModesR.Bot.rho;

% interpolate pressure onto the new grid
P = interp1( zL, PL, z, 'linear', 0.0 );

% special case for upper halfspace
if ( ModesR.Top.bc == 'A' )
    ii = find( z < depthTL );
    if ( ~isempty( ii ) )
        P( ii ) = phiTL(1:ML) * exp( -gamTL(1:ML) * ( depthTL - z( ii )' ) );
    end
end

% special case for lower halfspace
if ( ModesR.Bot.bc == 'A' )
    ii = find( z > depthBL );
    if ( ~isempty( ii ) )
        P( ii ) = phiBL(1:ML) * exp( -gamBL(1:ML) * ( z( ii )' - depthBL ) );
    end
end

% compute mesh width, h (it's actually h/rho)

% check this is correct for an irregular mesh

h = zeros( size( z ) );

h( 1           ) = 0.5 * ( z( 2       ) - z( 1           ) )  / rhomed( 1           ); % forward  diff
h( 2 : end - 1 ) = 0.5 * ( z( 3 : end ) - z( 1 : end - 2 ) ) ./ rhomed( 2 : end - 1 ); % centered diff
h( NTot        ) = 0.5 * ( z( NTot    ) - z( NTot - 1    ) )  / rhomed( NTot        ); % backward diff

% Point just above or below the interface
for med = 2 : ModesR.Nmedia
    % find points where the centered difference straddles an interface
    % typically 3 points will give two straddlers
    % hence the ./ rather than / below
    % see my notes
    DBelow  = ModesR.depth( med );
    ii = find( z( 1 : end - 2 ) < DBelow & z( 3 : end ) > DBelow ) + 1;
    
    h( ii ) = 0.5 * ( z( ii + 1 ) ./ rhomed( ii + 1 ) - z( ii - 1 ) ./ rhomed( ii - 1 ) ...
                         - DBelow ./ rhomed( ii + 1 ) +      DBelow ./ rhomed( ii - 1 ) );
end

P = h .* P;   % weight P by that h for later use in Trapezoidal rule


%**********************************************************************
function [ tail ] = tail( D, phiL, gamL, DL, phiR, gamR, DR )

% Computes the contributions from the tails in the halfspaces
% gamL and phiL are vectors
% tail = zeros( length( phiR ), 1 );

FL = phiL .* exp( -gamL * ( D - DL ) ).';
FR = phiR .* exp( -gamR * ( D - DR ) ).';

% gammaInv = 1 ./ ( gamL( :, ones( 1, length( gamR ) ) ) + gamR( :, ones( length( gamL ), 1 ) ).' );
gammaInv = 1 ./ ( gamL * ones( 1, length( gamR ) )  + ones( length( gamL ), 1 ) * gamR.' );

tail     = FR .* ( FL * gammaInv );
