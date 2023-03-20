function write_fieldsflp( flpfil, Pos )

% Write a field-parameters file

if ( ~strcmp( flpfil( end - 3 : end ), '.flp' ) )
  flpfil = [ flpfil '.flp' ]; % append extension
end

fid = fopen( flpfil, 'wt' );

rMin = min( Pos.r.r );
rMax = max( Pos.r.r );
NRr  = length( Pos.r.r );

fprintf( fid, '''RP'' \t \t ! Option \r\n ' );
fprintf( fid, '%6.2f %6.2f %5i \t ! rMin rMax (km) NRr \r\n', rMin, rMax, NRr );

fclose( fid );

