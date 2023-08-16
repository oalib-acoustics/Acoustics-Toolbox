% Runs a battery of test problems for the Acoustics Toolbox
% These are the ones out of the original KRAKEN manual
%
% note that KRAKEN will show lower TL result (higher levels) for the ICE case
% because it does not include attenuation in elastic media
%
% Similarly, SCOOTER shows a lower TL result for the kuperb case
% because it does not include attenuation due to scattering at interfaces. This includes the
% bottom interface

isd = 1;
ird = 1;

casesTL = [ 'pekeris' ; 'double '; 'scholte'; 'flused '; 'elsed  '; 'atten  '; ...
          'normal '; 'ice    '; 'kupert '; 'kuperb '; 'kuperm ' ];

copyfile( 'fieldbat.flp', 'field.flp' )
for icase = 1 : 11
   
   envfil = deblank( casesTL( icase, : ) );
   
   % kraken run
   copyfile( 'field.flp', [ envfil '.flp' ] );
   % disp( [ 'kraken( ''' envfil ''' )' ] );
   disp( envfil );

   eval( [ 'kraken( ''' envfil ''' )' ] );
   filename = [ envfil '.shd.mat' ];
   [ PlotTitle, PlotType, freq, freq0, atten, Pos, p ] = read_shd( filename );
   tlt( 1, : ) = -20.0 * log10( abs( p( ird, : ) ) );
   
   % krakenc run
   eval( [ 'krakenc( ''' envfil ''' )' ] );
   delete( [ envfil '.flp' ] )
   filename = [ envfil '.shd.mat' ];
   [ PlotTitle, PlotType, freq, freq0, atten, Pos, p ] = read_shd( filename );
   tlt( 2, : ) = -20.0 * log10( abs( p( ird, : ) ) );
   
   rkm = Pos.r.r / 1000;	% convert to km
   
   % plot
   
   figure
   plot( rkm, tlt( 1 : 2, : ) )
   xlabel( 'Range (km)' )
   ylabel( 'TL (dB)' );
   title( deblank( PlotTitle ) )
   
   set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
   hold on
   %plot( rkm, tlt( 1, : ), 'r-' )
   %plot( rkm, tlt( 2, : ), 'g-' )

   % scooter run
   if ( icase ~= 11 )   % scooter can't run kuperm with interfacial roughness
      copyfile( 'fields.flp', [ envfil '.flp' ] );
      eval( [ 'scooter( ''' envfil ''' )' ] );
      delete( [ envfil '.flp' ] )

      filename = [ envfil '.shd.mat' ]; % remove '.mat' to plot standard shdfil
      [ PlotTitle, PlotType, freq, freq0, atten, Pos, p ] = read_shd( filename );
      tlt( 3, : ) = -20.0 * log10( abs( p( ird, : ) ) );
      plot( Pos.r.r / 1000, tlt( 3, : ), 'k' )
      fclose( 'all' );
      title( deblank( PlotTitle ) )
   end
   
%    drawnow;
%    figure; plot( rkm, tlt( 1, : ) )
%    figure; plot( rkm, tlt( 2, : ) )
%    figure; plot( rkm, tlt( 3, : ) )

end

