function krakenM( filename )

% THIS IS NOT FINISHED !!!
% IT WILL NOT WORK

% Normal modes for ocean acoustics problems
% Based on the earlier Fortran version and Matlab version of Scooter
% MBP Sep. 12, 2006

clear global
tic

if ( isempty( filename ) )
    warndlg( 'No envfil has been selected', 'Warning' );
end

global omega Bdry
global NMedia N depth NFirstAcoustic NLastAcoustic H

% filenames
envfil = [ filename '.env' ];   % input environmental file
shdfil = [ filename '.mat' ];   % output file name (pressure)
trcfil = [ filename '.trc' ];   % input top reflection coefficient
brcfil = [ filename '.brc' ];   % input bottom reflection coefficient

[ PlotTitle, freq, Bdry, fid ] = readenv( envfil );   % read in the environmental file
PlotTitle = ['KRAKEN(M) -' PlotTitle ];

cLow   = fscanf( fid, '%f', 1 );   % lower phase speed limit
cHigh  = fscanf( fid, '%f', 1 );   % upper phase speed limit
fprintf( '\n cLow = %8.1f m/s   cHigh = %8.1f m/s  \n', cLow, cHigh )
fgetl( fid );

RMax  = fscanf( fid, '%f', 1 );   % read max range, Rmax
fprintf( 'RMax = %f \n', RMax )
fgetl( fid );

Pos = readsdrd( fid );            % read in the source and receiver depths
fclose( fid );                    % close out the envfil

readrc( trcfil, brcfil, Bdry.Top.Opt(2:2), Bdry.Bot.Opt(1:1) );   % READ Reflection Coefficients (top and bottom)

omega  = 2 * pi * freq;

% Set up vector of wavenumber samples
kMin = omega / cHigh;
kMax = omega / cLow;

NkPts = floor( 1000.0 * RMax * ( kMax - kMin ) / pi );
fprintf( '\nNumber of wavenumbers, NkPts = %i \n', NkPts )

H( 1:NMedia ) = ( depth( 2:NMedia + 1 ) - depth( 1:NMedia ) ) ./ N( 1:NMedia );   % vector of mesh widths
NptsAll       = sum( N( 1:NMedia ) ) + NMedia;   % number of solution points

[ B1, B2, B3, B4, rho ] = init_matrix( NptsAll );    % Initialize matrices
NptsAcoustic = sum( N( NFirstAcoustic:NLastAcoustic ) ) + 1;   % size of matrix for acoustic part
[ psi, k ] = init( B1, B2, B3, B4, rho, NptsAcoustic, Pos );
toc

save MODFIL PlotTitle freq psi k
%field( 'MODFIL.mat' );
%movefile( 'SHDFIL.mat', shdfil );
