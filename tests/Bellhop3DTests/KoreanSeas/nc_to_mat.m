% read the HYCOM oceanography files and convert to Matlab
% mbp 7/2011 based on earlier code

% either of the following commands can be used to show the variables in the
% netCDF file. ncdisp is built in to Matlab
%
% ncdump( 'ncom_relo_rimp_sml_2009080400.nc' )   % show the variables
% ncdisp archv.2011_204_00_3zt.nc

% Caution: the HYCOM grid is only regular below 40 degrees
% Therefore the part that extracts a subset will do odd things
% if your box goes above that

clear all

sspfil = 'Hycom3DSSP.asc';
files = [ ...
   'archv.2011_204_00_3' ];

% lat/lon box

latbox = [  30  41 ];
lonbox = [ 120 132 ];

% loop over analysis day

for iday = 1 : size( files, 1 )
   disp( iday )
   
   % open the temperature file
   ncid = netcdf.open( [ files( iday, : ) 'zt.nc' ], 'NC_NOWRITE' );
   
   % get the coordinates
   varid  = netcdf.inqVarID( ncid, 'Latitude' );   % latitude
   lat    = netcdf.getVar( ncid, varid );
   nlat   = size( lat, 2 );
   
   varid  = netcdf.inqVarID( ncid, 'Longitude' );   % longitude
   lon    = netcdf.getVar( ncid, varid );
   nlon   = size( lon, 1 );
   
   varid  = netcdf.inqVarID( ncid, 'Depth' );   % depth (m) positive down
   depth  = netcdf.getVar( ncid, varid );
   ndepth = length( depth );
   
   % get the temperature field
   varid = netcdf.inqVarID( ncid,        'temperature'   );
   t     = netcdf.getVar(   ncid, varid );
   
   % get the salinity file
   ncid  = netcdf.open( [ files( iday, : ) 'zs.nc' ], 'NC_NOWRITE' );
   varid = netcdf.inqVarID( ncid,        'salinity'   );
   s     = netcdf.getVar(   ncid, varid );

   % extract the user-selected lat/lon box
   
   jj = find( lat( 1, : ) >= latbox( 1 ) & lat( 1, : ) <= latbox( 2 ) );
   ii = find( lon( :, 2 ) >= lonbox( 1 ) & lon( :, 2 ) <= lonbox( 2 ) );
   
   lat = lat( 1, jj );
   lon = lon( ii, 1 );
   sbox = s( ii, jj );
   tbox = t( ii, jj );
   
   % tbox( tbox > 200 ) = 0; % NaN;   % replace fill values with NaNs
   % sbox( sbox > 100 ) = 0; % NaN;   % replace fill values with NaNs

   clear s t
   
   [ nx, ny ] = size( sbox )
   
   % convert to sound speed
   cmat = zeros( nx, ny, ndepth );
   for ii = 1 : nx
      for jj = 1 : ny
         cmat( ii, jj, : ) = soundspeed( squeeze( sbox( ii, jj, : ) ), squeeze( tbox( ii, jj, : ) ), depth );
      end
   end
   cmat = double( cmat );
   
   cmat( isnan( cmat ) ) = 1500.0;   % remove NaNs
   
   % convert the lat/lon coordinates to x-y (km)
   
   % Create grid in km
   % dist_wh returned a NaN for values on the same longitude (perturb by a
   % small number)
   [ LatTotkm, ~, A21 ] = dist_wh( [ lat( 1 ) lat( end ) ] , [ lon( 1 ) lon( 1   ) + .0001 ] );
   [ LonTotkm, ~, A21 ] = dist_wh( [ lat( 1 ) lat( 1   ) ] , [ lon( 1 ) lon( end ) ] );
   
   % mapping of lat/long to km
   Oceanography.x = linspace( 0, LonTotkm, length( lon ) ) / 1000;
   Oceanography.y = linspace( 0, LatTotkm, length( lat ) ) / 1000;
   
   Oceanography.lat   = lat;
   Oceanography.lon   = lon;
   Oceanography.depth = depth;
   %Oceanography.time  = tau;
   Oceanography.c( :, :, : ) = cmat;
   %%
   
   % plan view
   figure
   pcolor( Oceanography.lon, Oceanography.lat, squeeze( Oceanography.c( :, :, 5 ) )' ); shading flat;
   colorbar;
   caxis( [ 1500 1550 ] )
   xlabel( 'Longitude (degrees)' )
   ylabel( 'Latitude (degrees)' )
   
   figure
   pcolor( Oceanography.x, Oceanography.y, squeeze( Oceanography.c( :, :, 5 ) )' ); shading flat;
   colorbar;
   caxis( [ 1500 1550 ] )
   xlabel( 'x (km)' )
   ylabel( 'y (km)' )
   drawnow
   
   % side view
   lat_slice = 37;   % plot slice at this lattitude
   [ ~, jj ] = min( abs( Oceanography.lat - lat_slice ) );
   figure
   pcolor( Oceanography.x, Oceanography.depth, squeeze( Oceanography.c( :, jj, : ) )' );
   shading flat; colorbar;
   caxis( [ 1480 1540 ] )
   xlabel( 'x (km)' )
   ylabel( 'Depth (m)' )
   
   set( gca, 'YDir', 'Reverse' )   % plot with depth-axis positive down
   %%
   

   % save as a mat file
   save( [ files( iday, : ) '.mat' ], 'Oceanography' )
   
   % write to an ascii file for BELLHOP3D
   fid = fopen( sspfil, 'w' );
   fprintf( fid, '%i \r\n', nx );
   fprintf( fid, '%f  ', Oceanography.x );
   fprintf( fid, '\r\n' );
   
   fprintf( fid, '%i \r\n', ny );
   fprintf( fid, '%f  ', Oceanography.y );
   fprintf( fid, '\r\n' );
   
   fprintf( fid, '%i \r\n', ndepth );
   fprintf( fid, '%f  ', Oceanography.depth );
   fprintf( fid, '\r\n' );
   
   for kk = 1 : ndepth
      for jj = 1 : ny
         for ii = 1 : nx
            fprintf( fid, '%f  ', Oceanography.c( ii, jj, kk ) );
         end
         fprintf( fid, '\r\n' );
      end
   end
   fclose( fid );
   
end