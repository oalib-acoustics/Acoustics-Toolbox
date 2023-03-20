function krakelM( fileroot )

% Normal modes for ocean acoustics problems
% Based on the earlier Fortran version and Matlab version of Scooter
% MBP Sep. 12, 2006

clear global
tic

global Bdry SSP omega
global h

if ( isempty( fileroot ) )
    warndlg( 'No envfil has been selected', 'Warning' );
end

% filenames
envfil = [ fileroot '.env' ];       % input environmental file
%trcfil = [ fileroot '.trc' ];       % input top reflection coefficient
%brcfil = [ fileroot '.brc' ];       % input bottom reflection coefficient
modfil = [ fileroot '.mod.mat' ];   % output normal mode / eigenvalue, eigenvector file
%prtfil = [ fileroot '.prt' ];       % output print file

[ PlotTitle, freq, SSP, Bdry, ~, ~, ~, ~, fid ] = read_env( envfil, 'KRAKEL' );    % read in the environmental file
fclose( fid );                  % close out the envfil

% krakelM requires perfectly reflecting boundaries because it assumes an
% algebraic eigenvalue problem that can be passed to standard solvers
% You can emulate halfspaces with an absorbing layer of finite depth

if ( Bdry.Top.BC ~= 'V' && Bdry.Top.BC ~= 'R' )
   error( 'Fatal error: Top    Boundary Condition must be either V or R for krakelM' )
end

if ( Bdry.Bot.BC ~= 'V' && Bdry.Bot.BC ~= 'R' )
   error( 'Fatal error: Bottom Boundary Condition must be either V or R for krakelM' )
end

omega = 2 * pi * freq;

% Set up vector of wavenumber samples
% kMin = omega / cInt.High;
% kMax = omega / cInt.Low;

NMedia          = SSP.NMedia;
h( 1 : NMedia ) = ( SSP.depth( 2 : NMedia + 1 ) - SSP.depth( 1 : NMedia ) ) ./ SSP.N( 1 : NMedia );   % vector of mesh widths
NptsAll         = sum( SSP.N( 1 : NMedia ) ) + NMedia;   % number of mesh points

[ B1, B2, B3, B4, rho, Mater ] = init_elmatrix( NptsAll );   % initialize matrices

[ Asparse, Bsparse ] = setAB( B1, B2, B3, B4, rho, Mater );

% note: sparse routine 'eigs' requires Bsparse sym. p.d. It isn't.
%[ psi, k2 ] = eigs( Asparse, Bsparse, round( length( Asparse ) / 5 ), 'LR' );   %'LR' for largest real
[ psi, k2 ] = eig( full( Asparse ), full( Bsparse ) );

%Min_LOC = MINLOC( Extrap( 1, 1 : M ), Extrap( 1 : M ) > omega2 / CHigh^ 2 )
%M       = Min_LOC( 1 )

k = sqrt( diag( k2 ) );

% remove any eigenvalues with NaN's, infinities, or positive imaginary part
ii1  = find( ~isnan( k ) );
ii2  = find( ~isinf( k ) );
ii3  = find( imag( k ) <= 0 );

ii = intersect( ii1, ii2 );
ii = intersect( ii,  ii3 );

k   = k( ii );
psi = psi( :, ii );

% sort by real part of the eigenvalues
[ ~, ii ] = sort( real( k.^2 ), 1, 'descend' );
k   = k( ii );
psi = psi( :, ii );

% normalize each mode
M = length( k );   % number of modes

for Mode = 1 : M
  [ z, psi( :, Mode ) ] = normiz( psi( :, Mode ), k2( Mode ), B3, B4, rho, Mater );    % normalize
end

%psi( :, Mode ) = psi( 1 : length( z ), Mode );

% Note that the function normiz also calculates separate vectors u, v,
% tauzz, tauzz if you want those pulled out

% Write eigenvalues to PRTFil
fprintf( '\n\n     I          k             alpha          Phase Speed       Group Speed' )

for Mode = 1 : M
   fprintf( '\n %5i %16.10f %16.10f %16.10f', Mode, real( k( Mode ) ), imag( k( Mode ) ), omega / real( k( Mode ) ) )
end

fprintf( '\n\n' )
toc

% save the modes

Modes.title   = [ 'KRAKEL(M) - ', PlotTitle ];
Modes.Nfreq   = 1;
Modes.Nmedia  = SSP.NMedia;
Modes.N       = SSP.N;
Modes.depth   = SSP.depth;
Modes.Mater   = Mater;
Modes.rho     = rho;
Modes.freqVec = freq;    % freqVec( ifreq );
Modes.z       = z;
Modes.M       = M;
Modes.Top     = Bdry.Top;
Modes.Bot     = Bdry.Bot;
Modes.k       = k;
Modes.phi     = psi;

save( modfil, 'Modes' );
modfil
% calculate the pressure field
field( modfil );
