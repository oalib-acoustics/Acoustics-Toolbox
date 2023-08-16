function simplePE( filename )

% run the simplePE program
%
% usage: simplePE( filename )
% where filename is the environmental file (without the extension)

'NOTE!!! if the grid changes, then we need to update the SSP as well in RD cases'

% I have some mods in here that were just made to do a broadband
% calculation. So that is not in a very good state at the moment.
%
%
% (a simple MATLAB Parabolic equation model ...)
% M. Porter
% Uses a simple Crank-Nicholson scheme like the IFD PE
% See write up in JKPS
%
% This does *not* have 3 things that are useful in RAM:
% 1) does Pade (2,2) only, like the original IFD
% 2) does not do the Pade approximation to the full exponential operator
% ('split step)
% 3) does not do a partial LU decomposition that is re-usable
%
% The SimplePE *does* do linear interpolation of the SSP at each range step
%
% psi is the pressure (not reduced pressure)
%
% Complex sound speeds need to be conjugated relative to other models in the Acoustics
% Toolbox because of the different frequency convention
%
% 2021 August
% Fixed bug where sometimes the pressure field was not written at the last
% range.
% Fixed bug where the break statement was only exiting the receiver range
% loop and not the full PE range-marching loop
%
% First version written 7/94 at the Scripps Institution of Oceanography
% Extensively updated Feb. 2014

tic

%%
% Here are some additonal parameters you can set to control the PE run

% Think about the two-way travel through the bottom and the accumulated
% attenuation
% Typically 50 wavelengths is enough to kill off the return from the final
% bottom
Nlambda = 50; % number of wavelengths to extend the bottom depth to

% If you use fewer than 20 points/wavelength, you may see a bit of noise
Nsamples_per_wavelength = 10;   % finite difference sampling
%Nsamples_per_wavelength = 4;   % finite difference sampling (too coarse for at/tests/Gulf)
c0                      = 1500;   % reference sound speed
StarterType             = 'Greene';

%%
if ( isempty( filename ) )
    warndlg( 'No envfil has been selected', 'Warning' );
end

global omega Bdry
global Pos
global xBot

% filenames
envfil = [ filename '.env' ];   % input environmental file
shdfil = [ filename '.shd.mat' ];   % output file name (pressure)
btyfil = [ filename '.bty' ];   % input bathymetry
atifil = [ filename '.ati' ];   % input altimetry
sspfil = [ filename '.ssp' ];   % input ssp file (if range-dependent SSP)

% read the environmental file
[ PlotTitle, freq, SSP, Bdry, Pos, Beam, ~, ~, fid ] = read_env( envfil, 'BELLHOP' );
fclose( fid );                    % close out the envfil

Pos.r.r = 1000.0 * Pos.r.r; % convert km to m

TopATI = Bdry.Top.Opt( 5 : 5 );
readati( atifil, TopATI, Bdry.Top.depth, Beam.Box.r );    % read in the altimetry  file

BotBTY = Bdry.Bot.Opt( 2 : 2 );
readbty( btyfil, BotBTY, Bdry.Bot.depth, Beam.Box.r );    % read in the bathymetry file

is = 1;    % index of source depth
zs = Pos.s.z( is );			% source depth

freq_vec = linspace( 10, 100, 901 );
freq_vec = freq;   % this forces a run for a single frequency
nfreq    = length( freq_vec );

%%
% loop over frequencies

for ifreq = 1 : nfreq
    
    % if this is a broadband run, append the frequency index to the shdfil
    if ( nfreq > 1 )
        freq = freq_vec( ifreq );
        shdfil = [ filename int2str( ifreq ) '.shd.mat' ];   % output file name (pressure)
    end
    
    omega  = 2 * pi * freq;
    lambda = c0 / freq;
    k0     = omega / c0;
    
    %%
    % range grid
    % make sure the last receiver range is inside the grid
    rmax   = 1.01 * Pos.r.r( end );
    nr     = fix( Nsamples_per_wavelength * rmax / lambda );   % number of finite difference points
    deltar = rmax / nr;
    r      = linspace( deltar, rmax, nr );
    
    Nrd = length( Pos.r.z );
    Nrr = length( Pos.r.r );
    
    psit = zeros( Nrd, Nrr );   % sampled at receivers
    
    % get interpolated bottom depth at the r points of the f.d. grid
    d_interp = interp1( xBot( 1, : ), xBot( 2, : ), r );
    
    %%
    d = d_interp( 1 );   % bottom depth !!! assumes bathymetry starts at r=0 !!!
    d = max( d, lambda / 3 );   % make sure depth is not negative or ridiculously small

    % make initial depth grid
    % water
    nz = fix( Nsamples_per_wavelength * ( d - SSP.z( 1 ) ) / lambda );    % number of finite difference points
    nz = max( nz, 3 );

    h  = d / nz;  	% mesh spacing
    z  = SSP.z( 1 ) + linspace( h, d, nz );      % grid coordinates (surface point not part of the grid)
    
    % halfspace
    d_HS  = d + Nlambda * lambda;   % make absorbing layer Nlambda wavelengths thick
    nz_HS = fix( ( d_HS - d ) / h );   % number of finite difference points
    
    % adjust d_HS a bit so that the mesh size, h, divides it
    d_HS  = d + nz_HS * h;
    z_HS  = linspace( d + h, d_HS, nz_HS );   % f.d. grid for absorbing bottom layer
    
    z_tot  = [ z z_HS ];
    nz_tot = nz + nz_HS;
    
    %% optionally read a 2D SSP
    
    c = interp1( SSP.z, conj( SSP.c ), z ).';
    iProf = 1;   % counter to track the active profile at the current range step
    
    if ( Bdry.Top.Opt( 1 : 1 ) == 'Q' )
        [ cmat, rProf, NProf, ~ ] = readssp2d( sspfil );
        rProf = 1000 * rProf;   % convert km to meters
        % !!! assuming rProf( 1 ) = 0 here; should search for first profile
        c = interp1( SSP.z, cmat( :, iProf ), z ).';
    end
    
    %% Form marching matrix
    
    % generate ssp for absorbing bottom layer
    
    c_HS   = conj( Bdry.Bot.cp ) * ones( nz_HS, 1 );
    rho_HS = Bdry.Bot.rho;
    
    [ B, L, U ] = create_ABC( k0, c0, c, nz, d, c_HS, rho_HS, nz_HS, deltar );
    new_matrix = 0;   % flag to indicate whether the f.d. matrix needs updating
    
    
    %% March out in range
    
    ircvr_range = 1;
    
    psi = starter( StarterType, k0, zs, z_tot, nz_tot );
    
    for ir = 1 : nr
        range_march = r( ir );
        d_new       = d_interp( ir );
        d_new       = max( d_new, 0.25 );   % make sure depth is not negative or ridiculously small

        if ( mod( ir, fix( nr / 10 ) ) == 0 )
            fprintf( '      Range = %9.5g km \n', range_march / 1000 )
        end
        % if range-dependent bottom or SSP, update the matrices
        if ( BotBTY == '~'|| Bdry.Top.Opt( 1 : 1 ) == 'Q' )
            
            % show progress
            if ( mod( ir, fix( nr / 50 ) ) == 0 )
                fprintf( '      Range = %9.5g km   Depth = %9.5g m \n', range_march / 1000, d_new )
            end
            
            %% update the f.d. matrices if the depth changed or the SSP depends on range
            
            % NOTE!!!: the results are sensitive to the condition used in the following
            % if ( abs( d_new - d ) > .1 * lambda / Nsamples_per_wavelength )
            if ( abs( d_new - d ) > d / nz || Bdry.Top.Opt( 1 : 1 ) == 'Q' )
                [ z, nz, z_tot_new ] = make_grid( Nsamples_per_wavelength, d_new, lambda, nz_HS );
                
                % interpolate the pressure field onto the new grid
                % we add a point at the surface before doing interpolation

                psi   = interp1( [ 0 z_tot ], [ 0; psi ], z_tot_new, 'spline', 0 ).';
                z_tot = z_tot_new;
                d     = d_new;
                
                % c = interp1( SSP.z, SSP.c, z ).';
                c = interp1( SSP.z, conj( SSP.c ), z ).';

                new_matrix = 1;   % set flag to indicate we need to create a new marching matrix
            end
            
            %% update the f.d. matrix if the SSP changed or grid has changed
            
            if ( Bdry.Top.Opt( 1 : 1 ) == 'Q' )
                
                if ( range_march > rProf( iProf + 1 ) || mod( ir, fix( nr / 50 ) ) == 0  || new_matrix )
                    
                    % !!! assuming here there is a profile to find
                    ii = find( range_march < rProf );   % first profile left of current range
                    
                    if ( numel( ii ) > 0 )
                        iProf = ii( 1 ) - 1;
                        iProf = min( max( iProf, 1 ), NProf - 1 );
                        
                        c1 = interp1( SSP.z, cmat( :, iProf     ), z ).';
                        c2 = interp1( SSP.z, cmat( :, iProf + 1 ), z ).';
                        
                        % proportional distance in range between two SSP profiles
                        alpha = ( range_march - rProf( iProf ) ) / ( rProf( iProf + 1 ) - rProf( iProf ) );
                        
                        c     = ( 1 - alpha ) * c1 + alpha * c2;
                        new_matrix = 1;   % set flag to indicate we need to create a new marching matrix
                    end
                end
            end
            
            % Form marching matrix
            if ( new_matrix )
                [ B, L, U ] = create_ABC( k0, c0, c, nz, d, c_HS, rho_HS, nz_HS, deltar );
                new_matrix = 0;   % flag to indicate whether the f.d. matrix needs updating
            end
        end
        
        % equivalent to C psi_new = B * psi;
        psi_new = U \ ( L \ ( B * psi ) );
        
        %% interpolate the field onto the receiver grid
        % if we bracket a receiver, calculate the field at the receiver
        
        while range_march > Pos.r.r( ircvr_range )   % bracketted?
            
            deltar_temp = range_march - Pos.r.r( ircvr_range );  % range step to receiver
            
            if ( deltar_temp > 0 )
                % probably could avoid doing a new LU decomposition here
                
                [ Btemp, Ltemp, Utemp ] = create_ABC( k0, c0, c, nz, d, c_HS, rho_HS, nz_HS, deltar_temp );
                
                % equivalent to C psi_new = B * psi;
                psitemp = Utemp \ ( Ltemp \ ( Btemp * psi ) );
                
                psit( :, ircvr_range ) = interp1( z_tot, psitemp, Pos.r.z, 'linear', 0 );
                
                ircvr_range = ircvr_range + 1;
                if ( ircvr_range > Nrr )
                    break
                end   % jump out if marched past last receiver
            end
        end   % check next receiver range

        if ( ircvr_range > Nrr )
           break
        end   % jump out if marched past last receiver

        psi = psi_new;
    end   % next range step
    
    toc
    
    %% Scale and save the field
    
    % put back Hankel function
    rt = Pos.r.r';
    
    % scaled Hankel function matching (6.70) in JKPS
    hank = exp( 1i * ( k0 * rt - pi / 4 ) ) * diag( 1.0 ./ sqrt( rt ) );
    
    % conjugate the field to be consistent with the Acoustics Toolbox
    % convention
    is = 1;
    pressure( 1, is, :, : ) = conj( full( psit * spdiags( hank.', 0, Nrr, Nrr ) ) );
    
    PlotTitle    = ['SimplePE -' PlotTitle ];
    PlotType     = 'rectilin  ';
    atten        = 0;
    freqVec( 1 ) = freq;
    freq0        = freq;
    save( shdfil, 'PlotTitle', 'PlotType', 'freqVec', 'freq0', 'atten', 'Pos', 'pressure' )
    
end   % of frequency loop

end

%%
function [ B, L, U ] = create_ABC( k0, c0, c, nz, d, c_HS, rho_HS, nz_HS, deltar )

% sets up the matrices for the finite difference equations
% SPE: 2ik_0 { \pa \psi \over {\pa r} }
%          + { \pa^2 \psi \over {\pa z^2} }
%          + k_0^2 ( n^2 - 1 ) \psi = 0  (6.8)
%
% We have suppressed pivoting so that L and U are strictly tridiagonal

% ocean
n = c0 ./ c;
h = d / nz;	% mesh spacing

rho1 = 1;
rho2 = rho_HS;

% Absorbing bottom layer
n_HS   = c0 ./ c_HS;
nz_tot = nz + nz_HS;

% standard PE:
% E  = spdiags( ones( nz_tot, 1 ), -1, nz_tot, nz_tot );  % sub diagonal
% D1 = -2 * speye( nz_tot );  %     diagonal
% D2 = spdiags( [ h^2 * k0^2 * ( n.^2    - ones( nz,    1 ) );
%                 h^2 * k0^2 * ( n_HS.^2 - ones( nz_HS, 1 ) ) ], 0, nz_tot, nz_tot );
%
% A = D1 + D2 + E + E';
% B = 2 * 1i * k0 * h^2 / deltar * speye( size( A ) ) - A / 2;
% C = 2 * 1i * k0 * h^2 / deltar * speye( size( A ) ) + A / 2;
%

% Claerbout wide-angle PE:
a0 = 1.00;
a1 = 0.75;
b0 = 1.00;
b1 = 0.25;

w1  = b0 + .5 * 1i * k0 * deltar * ( a0 - b0 );
w1s = b0 - .5 * 1i * k0 * deltar * ( a0 - b0 );   % called w1* in jkps
w2  = b1 + .5 * 1i * k0 * deltar * ( a1 - b1 );
w2s = b1 - .5 * 1i * k0 * deltar * ( a1 - b1 );   % called w2* in jkps

r1 = w1  / w2;
r2 = w1s / w2s;

n2_tot = [ n.^2; n_HS.^2 ];

u    = 2 * ( k0^2 * h^2 / 2 * r2 - 1 ) + k0^2 * h^2 * ( n2_tot - 1 );
uhat = 2 * ( k0^2 * h^2 / 2 * r1 - 1 ) + k0^2 * h^2 * ( n2_tot - 1 );
nu = ones( nz_tot, 1 );

% jump at water/bottom interface
u( nz )    = ( rho1 + rho2 ) / rho2 * ( k0^2 * h^2 / 2 * r2 - 1 ) + ...
    .5 * k0^2 * h^2 * ( ( n2_tot( nz ) - 1 ) + rho1 / rho2 * ( n2_tot( nz + 1 ) - 1 ) );
uhat( nz ) = ( rho1 + rho2 ) / rho2 * ( k0^2 * h^2 / 2 * r1 - 1 ) + ...
    .5 * k0^2 * h^2 * ( ( n2_tot( nz ) - 1 ) + rho1 / rho2 * ( n2_tot( nz + 1 ) - 1 ) );
nu( nz + 1 )   = rho1 / rho2;   % offset by 1 to conform to Matlab spdiags convention

B = spdiags( ones( nz_tot, 1 ), -1, nz_tot, nz_tot ) + ...  % sub   diagonal
    spdiags( nu,                +1, nz_tot, nz_tot ) + ...  % super diagonal
    spdiags( uhat,               0, nz_tot, nz_tot );

% following equivalent is 2.5x slower
% B = gallery( 'tridiag', ones( nz_tot - 1, 1 ), uhat, nu( 2 : nz_tot ) );

B = w2 / w2s * B;

C = spdiags( ones( nz_tot, 1 ), -1, nz_tot, nz_tot ) + ...  % sub   diagonal
    spdiags( nu,                +1, nz_tot, nz_tot ) + ...  % super diagonal
    spdiags( u,                  0, nz_tot, nz_tot );

% C = gallery( 'tridiag', ones( nz_tot - 1, 1 ), u, nu( 2 : nz_tot ) );

thresh = 0.0;   % threshold for pivoting; 0 suppresses pivoting
[ L, U ] = lu( C, thresh );   % factor C

end

%% starter
function psi = starter( StarterType, k0, zs, z_tot, nz_tot )

switch ( StarterType )
    case( 'Gaussian' )
        % \psi(0,z) = \sqrt{ k_0 } \,
        % e^{ -{ k_0^2 \over 2 } ( z - z_s )^2 } (6.100)
        
        % Change fac to broaden the Gaussian, making it more narrow angle
        fac = 1;        % the usual formula has fac = 1
        
        psi = sqrt( k0 / fac ) * ...
            ( exp( -.5 * ( k0 / fac )^2 * ...
            ( ( z_tot - zs * ones( 1, nz_tot ) ) ).^2 )' ... % direct
            - exp( -.5 * ( k0 / fac )^2 * ...
            ( ( z_tot + zs * ones( 1, nz_tot ) ) ).^2 )' );  % surface image
        
    case( 'Greene' )
        % Greene's source JKPS(6.102)
        psi = sqrt( k0 ) * ( 1.4467 - 0.4201 * k0^2 * ( z_tot.' - zs * ones( nz_tot, 1 ) ).^2 ) .* ...
            ( exp( -( k0 )^2 * ( 1 / 3.0512 ) * ...
            ( ( z_tot.' - zs * ones( nz_tot, 1 ) ) ).^2 ) ) ... % direct
             -sqrt( k0 ) * ( 1.4467 - 0.4201 * k0^2 * ( z_tot.' + zs * ones( nz_tot, 1 ) ).^2 ) .* ...
            ( exp( -( k0 )^2 * ( 1 / 3.0512 ) * ...
            ( ( z_tot.' + zs * ones( nz_tot, 1 ) ) ).^2 ) ); % surface image  
         %figure; plot( z_tot, psi )
    case( 'delta' )
        % delta function starter:
        % psi = zeros( 1, nz_tot );
        %
        % [ ~, isd ] = min( z_tot - zs );
        % b   = zeros( 1, nz_tot );
        % b( isd ) = 1.;
    otherwise
        error( 'Unknown option for starting field' )
end
end

%%

function [ z, nz, z_tot ] = make_grid( Nsamples_per_wavelength, d, lambda, nz_HS )
% make the depth grid

% water
nz = fix( Nsamples_per_wavelength * d / lambda );    % number of finite difference points
nz = max( nz, 3 );

h  = d / nz;  	% mesh spacing
z  = linspace( h, d, nz );      % grid coordinates (surface point not part of the grid)

% adjust d_HS a bit so that the mesh size, h, divides it
d_HS  = d + nz_HS * h;
z_HS  = linspace( d + h, d_HS, nz_HS );   % f.d. grid for absorbing bottom layer

z_tot = [ z z_HS ];

% nz
% nz_HS
% nz_tot

end