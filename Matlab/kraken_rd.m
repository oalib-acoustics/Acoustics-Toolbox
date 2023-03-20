function kraken_rd( FileRoot, r )

% Calculate modes for a specified bathymetry
% This version is used with field.m (not fieldLoadAll.m)

% expects
% FileRoot.env as a Template
% FileRoot.bty containing the bathymetry
% r is a vector of ranges where the modes are calculated

global Bdry
global xBot

envfil = [ FileRoot '_rd' ];

% read in the template environment
model  = 'KRAKENC';
[ TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, ~ ] = read_env( FileRoot, model );

% optionally read bottom bathymetry

if ( Bdry.Bot.Opt(2:2) == '~' || Bdry.Bot.Opt(2:2) == '*' )
   btyfil = [ FileRoot '.bty' ];
   BotBTY = '~';
   depthB = SSP.depth( end );
   rBox   = RMax;
   readbty( btyfil, BotBTY, depthB, rBox )
end

D = interp1( xBot( 1, : ), xBot( 2, : ), r );

for ir = 1 : length( r )
   % substitute using the depth from the bathymetry file
   SSPtemp = SSP;   % make copy of template for temporary changes
   SSPtemp.depth( 2 ) = D( ir );
   
   % remove depths that exceed the bathymetry
   % need to be careful here in case D( ir ) is within rounding of the SSP
   % depth. The envfil only allows .01 m resolution ...
   ii = find( SSP.raw( 1 ).z > D( ir ) -.01, 1, 'first' );
   
   if ( ~isempty( ii ) )
      SSPtemp.raw( 1 ).z = SSP.raw( 1 ).z( 1 : ii - 1 );
      
      % add a new SSP point interpolated to the bathymetry
      if ( D( ir ) > SSP.raw( 1 ).z( ii - 1 ) )   % make sure added point is greater in depth
         SSPtemp.raw( 1 ).z(      ii ) = D( ir );
         SSPtemp.raw( 1 ).alphaR( ii ) = interp1( SSP.raw( 1 ).z, SSP.raw( 1 ).alphaR, D( ir ) );
      end
      
      % shift next layer up
      SSPtemp.raw( 2 ).z( 1   ) = D( ir );
   end
   
   TitleEnv2 = [ TitleEnv ' #' int2str( ir ) ];
   
   if ( ir == 1 )
      write_env( envfil, model, TitleEnv2, freq, SSPtemp, Bdry, Pos, Beam, cInt, RMax )
   else
      write_env( envfil, model, TitleEnv2, freq, SSPtemp, Bdry, Pos, Beam, cInt, RMax, 'a' )
   end
   % pause( ( length( r ) - ir ) / length( r ) * 1.0 )
end	% next range

% this is doing the same as kraken.m except that it allows me to force the
% use of the fortran field.f90

runkraken = which( 'kraken.exe' );

if ( isempty( runkraken ) )
   error( 'kraken.exe not found in your Matlab path' )
else
   eval( [ '! "' runkraken '" ' envfil ] );
end
