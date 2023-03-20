function bellhopM( filename )

% Beam-tracing code for ocean acoustics problems
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
% MBP Dec. 30, 2001
%
% For Matlab, the code had to be vectorized so all beams are traced in parallel.
% Furthermore, to keep the various vectors associated with each ray
% together in memory, it was necessary to collect them in a structure.
%
% The tracing of beams in parallel lead to huge storage requirements for the ray trajectories
% Suggests a possible rewrite to calculate TL (INFLUG module) in parallel
% with the ray trace (TRACE module).
% In the meantime, note that the run time goes through the ceiling when the
% required storage gets too large ...

clear global

tic

if ( isempty( filename ) )
    warndlg( 'No envfil has been selected', 'Warning' );
end

global ray MaxSteps Nsteps omega Bdry
global Pos
global alpha Nbeams
global U Arr MxNarr

% filenames
envfil = [ filename '.env' ];   % input environmental file
shdfil = [ filename '.shd.mat' ];   % output file name (pressure)
arrfil = [ filename '.arr.mat' ];   % output file name (arrivals)
btyfil = [ filename '.bty' ];   % input bathymetry
atifil = [ filename '.ati' ];   % input altimetry
% sbpfil = [ filename '.sbp' ];   % input source beam pattern
trcfil = [ filename        ];   % input top    reflection coefficient
brcfil = [ filename        ];   % input bottom reflection coefficient

[ PlotTitle, freq, SSP, Bdry, Pos, Beam, ~, ~, fid ] = read_env( envfil, 'BELLHOP' );    % read in the environmental file
fclose( fid );                    % close out the envfil

Beam.Box.r  = 1000.0 * Beam.Box.r;
Pos.r.r     = 1000.0 * Pos.r.r; % convert km to m

PlotTitle = ['BELLHOP(M) -' PlotTitle ];

SSP.Type = Bdry.Top.Opt(1:1);

% initialization for interpolation options other than c-Linear
switch ( SSP.Type )
   case ( 'N' )   % n^2 linear
      SSP.n2  = 1.0 ./ ( SSP.c.^2 );
      SSP.n2z = diff( SSP.n2 ) ./ diff( SSP.z );
   case ( 'P' )   % monotone PCHIP ACS (almost a cubic spline)
      SSP.ppc  = pchip_acs( SSP.z, SSP.c );
   case ( 'S' )   % Cubic Spline using not-a-knot boundary condition
      SSP.ppc = spline( SSP.z, SSP.c );
end

if ( SSP.Type == 'P'  || SSP.Type == 'S')
      [ breaks, coefs, pieces, order, dim ] = unmkpp( SSP.ppc );

      coefs_cz  = repmat( order - 1 : -1 : 1, pieces, 1 ) .* coefs( : , 1 : order - 1 );
      SSP.ppcz  = mkpp( breaks, coefs_cz,  dim ); % reassemble the differentiated piecewise polynomial

      coefs_czz = repmat( order - 2 : -1 : 1, pieces, 1 ) .* coefs_cz( : , 1 : order - 2 );
      SSP.ppczz = mkpp( breaks, coefs_czz, dim ); % reassemble the differentiated piecewise polynomial
end

TopATI = Bdry.Top.Opt(5:5);
readati( atifil, TopATI, Bdry.Top.depth, Beam.Box.r );    % read in the altimetry  file

BotBTY = Bdry.Bot.Opt(2:2);
readbty( btyfil, BotBTY, Bdry.Bot.depth, Beam.Box.r );    % read in the bathymetry file

readrc( trcfil, brcfil, Bdry.Top.Opt( 2 : 2 ), Bdry.Bot.Opt( 1 : 1 ) );   % READ Reflection Coefficients (top and bottom)
SrcBmPat = readpat( filename, Beam.RunType( 3 : 3 ) );   % optionally read source beam pattern

if ( ismember( 'A', upper( Beam.RunType( 1 : 1 ) ) ) )
    MxNarr = round( max( 2000000 / ( Pos.Nrz * Pos.Nrr ), 10 ) );   % allow space for at least 10 arrivals
    fprintf( '\n( Maximum # of arrivals = %i ) \n', MxNarr );
    
    Arr.A     = zeros( Pos.Nrz, Pos.Nrr, MxNarr );
    Arr.delay = zeros( Pos.Nrz, Pos.Nrr, MxNarr );
    Arr.angle = zeros( Pos.Nrz, Pos.Nrr, MxNarr );
    Arr.Narr  = zeros( Pos.Nrz, Pos.Nrr );
else
    pressure = zeros( 1, Pos.Nsz, Pos.Nrz, Pos.Nrr );
end

MaxSteps = round( 2000000 );

% pre-allocate
% ray( MaxSteps ) = struct( 'c', [], 'x', [], 'Tray', [ 0 0 ], 'p', [ 0 0 ], 'q', [ 0 0 ], 'tau', [], 'Rfa', [] );
ray( MaxSteps  ).c    = [];
ray( MaxSteps  ).x    = [];
ray( MaxSteps  ).Tray = [ 0 0 ];
ray( MaxSteps  ).p    = [ 0.0 0.0 ];
ray( MaxSteps  ).q    = [ 0.0 0.0 ];
ray( MaxSteps  ).tau  = [];
ray( MaxSteps  ).Rfa  = [];

omega = 2 * pi * freq;
Amp0  = interp1( SrcBmPat( :, 1 ), SrcBmPat( :, 2 ), alpha ); % calculate initial beam amplitude based on source beam pattern

alpha = ( pi / 180 ) * alpha;   % convert to radians
Dalpha = 0.0;
if ( Nbeams ~= 1 )
    Dalpha = ( alpha( Nbeams ) - alpha( 1 ) ) / ( Nbeams - 1 );  % angular spacing between beams
end

% *** Loop over source depths ***
Layer = 1;

for is = 1 : Pos.Nsz
    xs = [ 0.0 Pos.s.z( is ) ];   % source coordinate
    [ c, ~, ~, ~, ~, ~, Layer ] = ssp( xs, SSP, Layer );
    % RadMax = 10 * c / freq;  % 10 wavelength max radius
    
    % Are there enough beams?
    DalphaOpt = sqrt( c / ( 6.0 * freq * Pos.r.r( end ) ) );
    NbeamsOpt = round( 2 + ( alpha( Nbeams ) - alpha( 1 ) ) / DalphaOpt );
    
    if ( Beam.RunType(1:1) == 'C' && Nbeams < NbeamsOpt )
        fprintf( 'Warning: Nbeams should be at least = %i \n', NbeamsOpt )
    end
    
    % *** Trace rays ***
    U = zeros( length( Pos.r.z ), length( Pos.r.r ) );
    
    for ibeam = 1 : Nbeams
       
       if ( mod( ibeam - 1, max( fix( Nbeams / 50 ), 1 ) ) == 0 )
          disp( [' Tracing beam ', num2str( ibeam ) ] )
       end
       
        trace( SSP, Beam.deltas, xs, alpha( ibeam ), Amp0( ibeam ), Beam.Type, Beam.Box, Beam.RunType )
        
        % ray trace run?
        if Beam.RunType(1:1) == 'R'
            %figure
            title( PlotTitle );
            set( gca, 'YDir', 'Reverse' )   % plot with depth-axis positive down
            
            r = zeros( Nsteps, 1 );
            z = zeros( Nsteps, 1 );
            
            for ii = 1 : Nsteps
                r( ii ) = ray( ii ).x( 1 );
                z( ii ) = ray( ii ).x( 2 );
            end
            %if NumTopBnc >= 1 && NumBotBnc >= 1
            plot( r, z, 'k' )    % hits both boundaries
            %elseif NumBotBnc >= 1
            %    plot( r, z, 'b' )	   % hits bottom only
            %elseif NumTopBnc >= 1
            %    plot( r, z, 'g' )	   % hits surface only
            %else
            %    plot( r, z, 'r' )     % direct path
            %end
            drawnow; hold on
            %plot( r, z ); drawnow; hold on

            set( gca, 'YDir', 'Reverse' )   % plot with depth-axis positive down
            
            %save RAYFIL ibeams Nsteps ray
        else
            
            if ( Beam.RunType(2:2) == 'G' )
                InfluenceGeoHat(      Pos.s.z( is ), alpha( ibeam ), Beam.RunType, Dalpha  )
            else
                InfluenceGeoGaussian( Pos.s.z( is ), alpha( ibeam ), Beam.RunType, Dalpha, Beam.deltas  )
            end
        end
    end
    if ( Beam.RunType(1:1) == 'A' )   % arrivals calculation
        
        % add in cylindrical spreading
        for ir = 1 : Pos.Nrr
            factor = 1 / sqrt( Pos.r.r( ir ) );
            Arr.A( :, ir, : ) = factor * Arr.A( :, ir, : );
        end
        save( arrfil, 'PlotTitle', 'Pos', 'Arr', 'MxNarr' )
    else
        pressure( 1, is, :, : ) = scalep( Dalpha, c, Pos.r.r, U, Beam.RunType, Bdry.Top.Opt, freq );
        PlotType     = 'rectilin  ';
        atten        = 0;
        freqVec( 1 ) = freq;
        freq0        = freq;
        % problem: sd is outside the loop
        save( shdfil, 'PlotTitle', 'PlotType', 'freqVec', 'freq0', 'atten', 'Pos', 'pressure' )
        %plotshd( shdfil )
    end

end % next source depth

toc
