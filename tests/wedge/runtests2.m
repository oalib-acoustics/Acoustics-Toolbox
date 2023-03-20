% Run wedge problem to test coupled mode option in KRAKEN
global Bdry
global xBot tBot nBot tBotNode nBotNode RLenBot kappaBot NbtyPts btyType

Nprof = 51;
r = linspace( 0.0, 4000, Nprof );

% read in the template environment
model  = 'KRAKENC';
[ TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, fid ] = read_env( 'Template', model );

% optionally read bottom bathymetry
if ( Bdry.Bot.Opt(2:2) == '*' )
   btyfil = 'Template.bty';
   BotBTY = '*';
   depthB = SSP.depth( end );
   rBox   = RMax;
   readbty( btyfil, BotBTY, depthB, rBox )
end

D = interp1( xBot( 1, : ), xBot( 2, : ), r );

for ir = 1:length( r )
   envfil = 'ENVFIL';
   
   % substitute using the depth from the bathymetry file
   SSP.depth( 2 )      = D( ir );
   SSP.raw( 1 ).z( 2 ) = D( ir );
   SSP.raw( 2 ).z( 1 ) = D( ir );
   
   TitleEnv = [ 'Wedge problem #' int2str( ir ) ];

   if ( ir == 1 )
      write_env( envfil, model, TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax )
   else
      write_env( envfil, model, TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, 'a' )
   end
end	% next range

fclose( fid );

% run KRAKEN/FIELD

!krakenc.exe < ENVFIL > foo.prt
!field.exe < field.flp
%filename = [ 'Wedge' int2str( ifreq ) '.shd' ];

copyfile( 'SHDFIL', 'wedge.shd' );
delete MODFIL*

runplots
delete ENVFIL
