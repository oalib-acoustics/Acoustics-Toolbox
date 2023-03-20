function write_fieldflp( flpfil, Option, Pos )

% Write a field-parameters file

if ( ~strcmp( flpfil( end - 3 : end ), '.flp' ) )
  flpfil = [ flpfil '.flp' ]; % append extension
end

fid = fopen( flpfil, 'wt' );

fprintf( fid, '/ ! Title \n' );
fprintf( fid, '''%-4s''  ! Option \n', Option );
fprintf( fid, '999999   ! Mlimit (number of modes to include) \n' );
fprintf( fid, '1        ! NProf  \n' );
fprintf( fid, '0.0 /    ! rProf (km) \n' );

% receiver ranges
fprintf( fid, '%5i \t \t \t \t ! NRr \n', length( Pos.r.r ) );

if ( length( Pos.r.r ) > 2 && equally_spaced( Pos.r.r ) )
    fprintf( fid, '    %6f  ', Pos.r.r( 1 ), Pos.r.r( end ) );
else
    fprintf( fid, '    %6f  ', Pos.r.r );
end
fprintf( fid, '/ \t ! Rr(1)  ... (km) \n' );

% source depths

fprintf( fid, '%5i \t \t \t \t ! NSz \n', length( Pos.s.z ) );

if ( length( Pos.s.z ) > 2 && equally_spaced( Pos.s.z ) )
    fprintf( fid, '    %6f  ', Pos.s.z( 1 ), Pos.s.z( end ) );
else
    fprintf( fid, '    %6f  ', Pos.s.z );
end

fprintf( fid, '/ \t ! Sz(1)  ... (m) \n' );

% receiver depths

fprintf( fid, '%5i \t \t \t \t ! NRz \n', length( Pos.r.z ) );

if ( length( Pos.r.z ) > 2 && equally_spaced( Pos.r.z ) )
    fprintf( fid, '    %6f  ', Pos.r.z( 1 ), Pos.r.z( end ) );
else
    fprintf( fid, '    %6f  ', Pos.r.z );
end

fprintf( fid, '/ \t ! Rz(1)  ... (m) \n' );

% receiver range offsets

fprintf( fid, '%5i \t \t \t \t ! NRro \n', length( Pos.r.z ) );
fprintf( fid, '    %6.2f  ', zeros( 1, 2 ) );
fprintf( fid, '/ \t \t \t \t ! Rro(1)  ... (m) \n' );

fclose( fid );
