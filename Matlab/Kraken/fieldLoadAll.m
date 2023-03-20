function fieldLoadAll( fileroot, Dwedge )

% calculates the field using modes produced by KRAKEN
% parameters read from field.flp
%
% usage: fieldLoadAll( fileroot )
% mbp

filename = [ fileroot '_rd.mod' ];
shdfil = [ fileroot '.shd.mat' ];

[ TitleEnv, Opt, Comp, MLimit, NProf, rProf, R, Pos ] = read_flp( fileroot );
modes = 1 : MLimit;
rd    = Pos.r.z;
Nrd   = length( rd );
Nr    = length( R );

clear read_modes_bin % to force rewind to beginning of mode file

% optionally read bottom bathymetry
%if ( Bdry.Bot.Opt(2:2) == '*' )
btyfil = [ fileroot '.bty' ];
BotBTY = '~';
depthB = 3000; % SSP.depth( end );
rBox   = 30; % RMax;
readbty( btyfil, BotBTY, depthB, rBox )
%end

if ( NProf == 1 )   % Range-independent case
   
   if nargin == 1
      Modes = read_modes( filename );
   else
      Modes = read_modes( filename, modes );
   end
   MSrc = length( Modes.k );
   M   = min( MLimit, MSrc );   % Set number of propagating modes

   % weight modes by mode excitation
   zs  = Pos.s.z;
   isd = find( Modes.z >= zs );		% index of source depth
   isd = isd( 1 );
   
   C   = Modes.phi( isd, 1 : min( MLimit, Modes.M ) ).';   % column vector with modal weights
   
   % calculate mode values at receiver depths
   irdvec = zeros( 1, length( Pos.r.z ) );
   for ii = 1 : length( Pos.r.z )
      zr  = Pos.r.z( ii );
      ird = find( Modes.z >= zr );		% index of source depth
      irdvec( ii ) = ird( 1 );
   end
   phiR = Modes.phi( irdvec, : );
   M   = min( MLimit, MSrc );   % Set number of propagating modes
   
   pressure = evalri( C, phiR( :, 1 : min( MLimit, Modes.M ) ), R, Modes.k( 1 : min( MLimit, Modes.M ) ), Opt, Comp );
else
   if ( Opt(2:2) == 'C' )   % Range-dependent case
      % Coupled mode theory
      pressure = evalcmLoadAll( fileroot, Pos, Dwedge, rd, Nrd, R, Nr, MLimit, Opt );
   else
      % Adiabatic mode theory
      pressure = evalad( rProf, NProf,      C, phiR, rd, Nrd, R, Nr, M, Opt );
   end
end

PlotTitle = TitleEnv;
PlotType  = 'rectilin  ';
freqVec   = 0; % Modes.freq;
atten     = 0.0;
Pos.r.r   = R;

save( shdfil, 'PlotTitle', 'PlotType', 'freqVec', 'atten', 'Pos', 'pressure' )

