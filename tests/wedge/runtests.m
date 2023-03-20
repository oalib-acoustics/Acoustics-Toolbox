function runtests
% Run wedge problem to test coupled mode option in KRAKEN
% write out the envfil

envfil   = 'wedge';
freq  = 100:1:200;
freq  = 25;
nfreq = length( freq );

Nprof    = 51;
r        = linspace( 0.0, 4000, Nprof );
D        = linspace( 200,    0, Nprof );
D( end ) = 1.0;

for ifreq = 1 : nfreq
   
   for ir = 1 : length( r )
      fprintf( ' Range index %4i    Range %7.1f \n', ir, r( ir ) )
      model    = 'KRAKENC';
      TitleEnv = [ 'Wedge problem #' int2str( ir ) ];
      freqT    = freq( ifreq );
      
      SSP.NMedia = 2;
      SSP.N      = [ 500 1000 ];
      SSP.sigma  = [   0    0 ];
      SSP.depth  = [ 0 D( ir ) 2000 ];
      
      SSP.raw( 1 ).z      = [ 0 D( ir ) ];
      SSP.raw( 1 ).alphaR = [ 1500 1500 ];
      SSP.raw( 1 ).betaR  = [ 0 0 ];
      SSP.raw( 1 ).rho    = [ 1 1 ];
      SSP.raw( 1 ).alphaI = [ 0 0 ];
      SSP.raw( 1 ).betaI  = [ 0 0 ];
      
      SSP.raw( 2 ).z      = [ D( ir ) 2000 ];
      SSP.raw( 2 ).alphaR = [ 1700 1700 ];
      SSP.raw( 2 ).betaR  = [ 0   0   ];
      SSP.raw( 2 ).rho    = [ 1.5 1.5 ];
      SSP.raw( 2 ).alphaI = [ 0.5 0.5 ];
      SSP.raw( 2 ).betaI  = [ 0   0 ];
      
      Bdry.Top.Opt   = 'CVW .';
      Bdry.Bot.Opt   = 'V';
      Bdry.Bot.sigma = 0.0;
      
      cInt.Low  = 1400;
      cInt.High = 15000;
      
      RMax    = 0.0;
      Pos.s.z = 0.0;
      Pos.r.z = linspace( 0, 2000, 2001 );
      
      Beam = 'dummy';
      
      if ( ir == 1 )
         write_env( envfil, model, TitleEnv, freqT, SSP, Bdry, Pos, Beam, cInt, RMax )
      else
         write_env( envfil, model, TitleEnv, freqT, SSP, Bdry, Pos, Beam, cInt, RMax, 'a' )
      end
   end	% next range
      
   % run KRAKEN/FIELD
   
   runkraken = which( 'krakenc.exe' );
   eval( [ '! "' runkraken '" ' envfil ] );
   
   runfield = which( 'field.exe' );
   eval( [ '! "' runfield '" ' envfil ] );

end   % next frequency

runplots

%%

% wedge flat
% what to expect?
% steep angles will be off because of the narrow angle approximation
% a little noise from the bottom reflection

scooter wedge_flatS
plotshd( 'wedge_flatS.shd.mat', 2, 1, 1 )
caxisrev( [ 30 80 ] )

simplePE wedge_flatB
plotshd( 'wedge_flatB.shd.mat', 2, 1, 2 )
caxisrev( [ 30 80 ] )

% wedge
% agreement with coupled modes should be nearly perfect
simplePE wedgeB
figure
plotshd( 'wedgeB.shd.mat' )
caxisrev( [ 40 80 ] )

