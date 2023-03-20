function writebdry( bdryfil, interp_type, rngdep )

% Write a boundary file (bathymetry or altimetry) from the workspace variables
% ranges should be in km
% depths should be in m

Npts = length( rngdep( :, 1 ) );

fid = fopen( bdryfil, 'wt' );
fprintf( fid, '''%c''', interp_type );
fprintf( fid, '\n');

fprintf( fid, '%i', Npts );
fprintf( fid, '\n');

fprintf( fid, '%f %f \n', rngdep' );
fprintf( fid, '\n');

fclose( fid );

