function sparcM( filename )

% Wavenumber integration code for ocean acoustics problems
% Copyright (C) 2009 Michael B. Porter
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Based on the earlier Fortran version
% MBP Dec. 24, 2005

% Note: !!! This Matlab version does not have the source filtering that is
% done in the Fortran version. The source filtering forms a bandpass
% waveform, pre-envelope, and/or Hilbert transform of s(t).
% This allows the Fortran version, for instance, to remove the
% left-traveling wave.

clear global
tic

if ( isempty( filename ) )
    warndlg( 'No envfil has been selected', 'Warning' );
end

global omega Bdry
global H
global Pulse fMin fMax tStart tMult alpha beta V
global tout

% filenames
envfil = [ filename '.env' ];   % input environmental file

[ PlotTitle, freq, SSP, Bdry, Pos, ~, cInt, RMax, fid ] = read_env( envfil, 'SPARC' );    % read in the environmental file

PlotTitle = ['SPARC(M) -' PlotTitle ];
omega  = 2 * pi * freq;

% *** souce information ***

Pulse = fgetl( fid );

% Extract option letters between the quotes
nchars = strfind( Pulse, '''' );   % find quotes
Pulse   = [ Pulse( nchars( 1 ) + 1 : nchars( 2 ) - 1 ) blanks( 5 - ( nchars( 2 ) - nchars( 1 ) ) ) ];

switch ( Pulse(1:1) )
    case ( 'P' )
        fprintf( 'Pseudo-gaussian pulse' )
    case ( 'R' )
        fprintf( 'Ricker wavelet' )
    case ( 'A' )
        fprintf( 'Approximate Ricker wavelet' )
    case ( 'S' )
        fprintf( 'Single sine source' )
    case ( 'H' )
        fprintf( 'Hanning weighted four sine pulse' )
    case ( 'N' )
        fprintf( 'N-wave pulse' )
    case ( 'M' )
        fprintf( 'Miracle-wave pulse' )
    case ( 'G' )
        fprintf( 'Gaussian pulse' )
    case ( 'F' )
        fprintf( 'Source time series from File' )
    case ( 'B' )
        fprintf( 'Source time series reversed from file' )
    otherwise
        error( 'Unknown source type' )
end

fMin  = fscanf( fid, '%f', 1 );
fMax  = fscanf( fid, '%f', 1 ); % Upper and lower frequency limits
fgetl( fid );
fprintf( '\nfMin, fMax, %f %f', fMin, fMax )

kMin = 2.0 * pi * fMin / cInt.High;   % minimum wavenumber
kMax = 2.0 * pi * fMax / cInt.Low;    % maximum wavenumber
if ( cInt.High > 1.0E6 )
    kMin = 0.0;
end

Nk = 1000.0 * RMax * ( kMax - kMin ) / ( 2.0 * pi );
fprintf( '\n\nNumber of wavenumbers, Nk = %i \n', Nk )
k = linspace( kMin, kMax, Nk );

rr = readr( fid ); % Read receiver ranges
rr = 1000.0 * rr; % convert km to m

tout = readt( fid ); % Read in the output times

% *** Integration parameters ***
% alpha = 0.0  lumped     mass matrix,
%         1.0  consistent mass matrix
% beta  = 0.0  standard explicit,
%         0.25 standard implicit

tStart = fscanf( fid, '%f', 1 );
tMult  = fscanf( fid, '%f', 1 );
alpha  = fscanf( fid, '%f', 1 );
beta   = fscanf( fid, '%f', 1 );
V      = fscanf( fid, '%f', 1 );

fprintf( '\n tStart = %f', tStart );
fprintf( '\n tMult  = %f', tMult );
fprintf( '\n alpha  = %f', alpha );
fprintf( '\n beta   = %f', beta );
fprintf( '\n V      = %f', V );
fclose( fid );      % close out the envfil

H( 1:SSP.NMedia ) = ( SSP.depth( 2:SSP.NMedia + 1 ) - SSP.depth( 1:SSP.NMedia ) ) ./ SSP.N( 1:SSP.NMedia );   % vector of mesh widths
[ c2R, c2I, rho, crosst, cMin, cMax ] = initsparc( SSP, Pos );

kernelsparc( filename, SSP, c2R, c2I, rho, crosst, cMin, cMax, k, freq, Pos, rr, Bdry, PlotTitle );

fprintf( '\n' );
toc