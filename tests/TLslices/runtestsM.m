% Runs a battery of test problems for the Acoustics Toolbox
% These are the ones out of the original KRAKEN manual

isd = 1;
ird = 1;

cases = [ 'pekeris' ; 'double '; 'scholte'; 'flused '; 'elsed  '; 'atten  '; 'normal '; 'ice    '; 'kupert '; 'kuperb '; 'kuperm ' ] 

copyfile( 'fieldbat.flp', 'field.flp' )
for icase = 1:11
   
   envfil = deblank( cases( icase, : ) )
   
   % kraken run
   
   eval( [ 'kraken( ''' envfil ''')' ] );
   filename = [ envfil '.shd' ];
   [ PlotTitle, PlotType, freq, atten, Pos, p ] = read_shd( filename );
   tlt( 1, : ) = -20.0 * log10( abs( p( ird, : ) ) );
   
   % krakenc run
   eval( [ 'krakenc( ''' envfil ''' )' ] );
   filename = [ envfil '.shd' ];
   [ PlotTitle, PlotType, freq, atten, Pos, p ] = read_shd( filename );
   tlt( 2, : ) = -20.0 * log10( abs( p( ird, : ) ) );
   
   rkm = Pos.r.range / 1000;	% convert to km
   % plot
   
   figure
   plot( rkm, tlt(1:2, : ) )
   xlabel( 'Range (km)' ); ylabel( 'TL (dB)' );
   title( deblank( PlotTitle ) )
   set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
   hold on
   
   % scooter run
   if ( icase ~= 11 )   % scooter can't run kuperm with interfacial roughness
   eval( [ 'scooterM( ''' envfil ''' )' ] );
   filename = [ envfil '.mat' ]; % remove '.mat' to plot standard shdfil
   [ PlotTitle, PlotType, freq, atten, Pos, p ] = read_shd( filename );
   tlt( 3, : ) = -20.0 * log10( abs( p( ird, : ) ) );
   plot( Pos.r.range/1000, tlt( 3, : ) )
   fclose( 'all' )
   title( deblank( PlotTitle ) )
   end
   
   drawnow;
   
end

