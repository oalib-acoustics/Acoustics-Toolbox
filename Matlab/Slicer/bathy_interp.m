function rd = bathy_interp( Bathy, Interp )

% BATHY_INTERP - interpolate a bathymetry struct
% bearings are in degrees using the math convention

% extract the needed fields from the Interp struct

src_lat = Interp.src_lat;
src_lon = Interp.src_lon;

bearing_deg    = Interp.bearing_deg;
bathy_range_km = Interp.bathy_range_km;
bathy_delta_km = Interp.bathy_delta_km;

% compute the number of range, depth points and tweak max range

npt = 1 + ceil( bathy_range_km / bathy_delta_km );

bathy_range_km = ( npt - 1 ) * bathy_delta_km;

% compute the (Northing, Easting) values for the start of the plane

beg_northing_km = interp1( Bathy.Lat, Bathy.Latkm, src_lat );
beg_easting_km  = interp1( Bathy.Lon, Bathy.Lonkm, src_lon );

% compute the (Northing, Easting) values for the endpoint of the plane

bearing_rad = ( pi / 180.0 ) * bearing_deg;

end_northing_km = beg_northing_km + bathy_range_km * sin( bearing_rad );
end_easting_km  = beg_easting_km  + bathy_range_km * cos( bearing_rad );

% interpolate bathymetry at equally spaced range points

range_km    = linspace( 0.0,             bathy_range_km,  npt );
easting_km  = linspace( beg_easting_km,  end_easting_km,  npt );
northing_km = linspace( beg_northing_km, end_northing_km, npt );

depth_m = interp2( Bathy.Lonkm, Bathy.Latkm, Bathy.depth, ...
                   easting_km, northing_km, 'linear', 0.0 );

% arrange the interpolated range and depth values into an Nx2 matrix

rd = [ range_km( : ), depth_m( : ) ];
