function write_bell( fid, Beam )

% Write out the rest of the environmental file
% This is the part containing control info specific to Bellhop

% *** write run type ***
fprintf( fid, '''%s'' \t \t \t \t ! Run Type \n', Beam.RunType );

if ( isfield( Beam, 'Ibeam' ) )
   fprintf( fid, '%i %i \t \t \t \t ! Nbeams Ibeam \n', Beam.Nbeams, Beam.Ibeam );
else
   % if this is a ray trace run and the field Beam.Nrays exists to use
   % fewer rays in the trace, then use that
   if ( Beam.RunType( 1 : 1 ) == 'R' && isfield( Beam, 'Nrays' ) )
      fprintf( fid, '%i \t \t \t \t \t! Nbeams \n', Beam.Nrays );
   else
      fprintf( fid, '%i \t \t \t \t \t! Nbeams \n', Beam.Nbeams );
   end
end

fprintf( fid, '%f %f / \t \t ! angles (degrees) \n', Beam.alpha( 1 ), Beam.alpha( end ) );
fprintf( fid, '%f %f %f \t ! deltas (m) Box.z (m) Box.r (km) \n', Beam.deltas, Beam.Box.z, Beam.Box.r );

% Cerveny-style Gaussian beams
if ( length( Beam.RunType ) > 1 &&  ~contains( 'GBS', Beam.RunType(2:2) )  )
   fprintf( fid, '''%s'' %f %f  \t \t ! ''Min/Fill/Cer, Sin/Doub/Zero'' Epsmult RLoop (km) \n', Beam.Type(1:2), Beam.epmult, Beam.rLoop );
   fprintf( fid, '%i %i  \t \t \t \t ! Nimage Ibwin \n', Beam.Nimage, Beam.Ibwin );
end