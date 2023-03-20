function scooterM( filename )

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

% type:
% mex factortri.c
% mex backsub.c
% to have Matlab call these c-versions in place of factortri.m and
% backsub.m
% The speed up is 5x in the last test I did
% The mex files are not generally portable between different architectures,
% hence the need to do the mex command yourself

% Based on the earlier Fortran version
% MBP Dec. 24, 2003

% clear global
tic
global omega SSP Bdry
global NFirstAcoustic NLastAcoustic H

if ( isempty( filename ) )
    warndlg( 'No envfil has been selected', 'Warning' );
end

% filenames
envfil = [ filename '.env' ];     % input environmental file
grnfil = [ filename '.grn.mat' ];
trcfil = filename;                % input top reflection coefficient
brcfil = filename;                % input bottom reflection coefficient

[ PlotTitle, freq, SSP, Bdry, Pos, ~, cInt, RMax, fid ] = read_env( envfil, 'SCOOTER' ); % read in the environmental file
fclose( fid );                    % close out the envfil

PlotTitle = ['SCOOTER(M) -' PlotTitle ];

readrc( trcfil, brcfil, Bdry.Top.Opt(2:2), Bdry.Bot.Opt(1:1) );   % READ Reflection Coefficients (top and bottom)

omega  = 2 * pi * freq;

% Set up vector of wavenumber samples
kMin = omega / cInt.High;
kMax = omega / cInt.Low;

NkPts = floor( 1000.0 * RMax * ( kMax - kMin ) / pi );
fprintf( '\nNumber of wavenumbers, NkPts = %i \n', NkPts )

% Set-up the vector of k-space points
DeltaK = ( kMax - kMin ) / ( NkPts - 1 );
atten  = DeltaK;
rk     = linspace( kMin, kMax, NkPts );

H( 1 : SSP.NMedia ) = ( SSP.depth( 2 : SSP.NMedia + 1 ) - SSP.depth( 1 : SSP.NMedia ) ) ./ SSP.N( 1 : SSP.NMedia ); % vector of mesh widths
NptsAll       = sum( SSP.N( 1 : SSP.NMedia ) ) + SSP.NMedia;   % number of solution points

[ B1, B2, B3, B4, rho ] = init_matrix( NptsAll );    % Initialize matrices

NptsAcoustic = sum( SSP.N( NFirstAcoustic : NLastAcoustic ) ) + 1;   % size of matrix for acoustic part
Green = kernel( B1, B2, B3, B4, rho, rk, atten, NptsAcoustic, Pos );
toc

PlotType = 'rectilin  ';
Pos.Nrr  = NkPts;
Pos.r.r  = omega ./ rk.';   % store phasespeeds in slot usually used for receiver ranges

% note pressure is 4-d: pressure( Ntheta, Nsd, Nrd, Nrr )
pressure( 1, 1 : Pos.Nsz, 1 : Pos.Nrz, 1 : Pos.Nrr ) = Green;  % store Green's function in slot usually used for pressure

freqVec( 1 ) = freq;
freq0        = freq;
save( grnfil, 'PlotTitle', 'PlotType', 'freqVec', 'freq0', 'atten', 'Pos', 'pressure' )
fieldsco( grnfil );
