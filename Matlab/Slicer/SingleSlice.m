function [ Slice_out, rv, rc ] = SingleSlice( ConfigParams, Slice, Bathy )
% COMPUTE_SLICE - 
%
% setup the information needed for Arrivals, or TL calculation and
% run the selected acoustics model for the given slice (vertical plane)
%
%    rv is the returned value (usually TL)
%    rc is a completion code (0 for success)
%    Slice_out is a modified version of Slice (the only thing that matters is
%       Slice_out.Pos.r.r, used by TL_grid)
%
% The primary tasks are to extract the bathymetry and range-dependent
% sound speed profile for the specified vertical plane.
%
% There are two methods available to specify the vertical plane. Provide:
%  1) lat/lon of the source and the lat/lon of most distant receiver
%  2) lat/lon of the source and the bearing and range scalar values

Slice_out = Slice;   % copy the input Slice struct (since it will be modified)

% extract the needed fields from the input structures

Plane           = Slice.Plane;
bathy_delta_km  = Slice.bathy_delta_km;
bty_file        = Slice.Bdry.Bot.Opt( 2 );

if isfield( Plane, 'min_d_m' )
   min_d_m = Plane.min_d_m;
else
   min_d_m = 0.0;
end

max_d_m   = Plane.max_d_m;
nrr       = Plane.nrr;
nrd       = Plane.nrd;
ssp_depth = Slice.SSP.depth( 2 );

% Compute the parameters that describe the computational plane

switch upper( Plane.method )
   case 'ENDPOINTS'
      src_lat = Plane.src_lat;
      src_lon = Plane.src_lon;
      
      rcv_lat = Plane.rcv_lat;
      rcv_lon = Plane.rcv_lon;
      
      % Compute the bearing and range of the receiver (relative to src)
      
      [ r_m, bearing_deg, ~ ] = dist_wh_math( [ src_lat rcv_lat ], [ src_lon rcv_lon ] );

      min_r_km = 0.0;
      max_r_km = r_m / 1000.0;
      
   case 'POLAR'
      src_lat = Plane.src_lat;
      src_lon = Plane.src_lon;
      
      bearing_deg = Plane.bearing_deg;
      
      if isfield( Plane, 'min_r_km' )
         min_r_km = Plane.min_r_km;
      else
         min_r_km = 0.0;
      end
      
      max_r_km = Plane.max_r_km;
      
   otherwise
      error( [ mfilename, ': unrecognized method for specification of computational plane: ', method ] );
end

% Setup the receiver positions on the computational plane

del_r = max_r_km / nrr;

if nrr == 1 || min_r_km < del_r
   rr = ( 1 : nrr ) * del_r;
else
   rr = linspace( min_r_km, max_r_km, nrr );
end

del_z = max_d_m  / nrd;

if nrd == 1 || min_d_m < del_z
   rd = ( 1 : nrd ) * del_z;
else
   rd = linspace( min_d_m, max_d_m, nrd );
end

Pos           = Slice.Pos;
Pos.Nrr       = nrr;
Pos.Nrd       = nrd;
Pos.r.r   = rr;
Pos.r.z   = rd;
Slice_out.Pos = Pos;

Slice_out.Beam.Box.r = max_r_km;

% If requested, interpolate the bathymetry along the computational plane

if bty_file == '*'
   
   % Determine range limit for bathymetry interpolation
   
   if nrr > 1
      bathy_range_km = rr( end ) + del_r;
   else
      bathy_range_km = rr( 1   ) + 0.5;
   end
   
   % Interpolate the bathymetry for the specified computational plane
   
   Interp.src_lat        = src_lat;
   Interp.src_lon        = src_lon;
   Interp.bearing_deg    = bearing_deg;
   Interp.bathy_range_km = bathy_range_km;
   Interp.bathy_delta_km = bathy_delta_km;
   
   rd = bathy_interp( Bathy, Interp );
   
   Slice_out.bathy_rd = rd;
   
   % Verify the SSP extends to the deepest point in the interpolated bathymetry
   
   max_depth = max( rd( :, 2 ) );
   
   if max_depth > ssp_depth
      error_msg = sprintf( ': max depth of sound speed profile (%.2f) does not extend to bottom (%.2f)', ...
         ssp_depth, max_depth );
      error( [ mfilename, error_msg ] );
   end
   
   % Reset the limits of the bounding box for ray tracing models (BELLHOP)
   
   Slice_out.Beam.Box.r = bathy_range_km;
   Slice_out.Beam.Box.z = max_depth + 1;
   
end

%%

% run the acoustic model

acoustic_model = upper( Slice_out.acoustic_model );

switch acoustic_model
   case 'BELLHOP'
      write_bellhop_files( ConfigParams, Slice_out );  % write the BELLHOP input files
      [ rv, rc ] = run_bellhop( ConfigParams, Slice_out );
      
   case 'RAM'
      write_ram_files( ConfigParams, Slice_out );  % write the RAM input files
      [ rv, rc ] = run_ram( ConfigParams );
      
   otherwise
      error( [ mfilename, ': error: unknown or unsupported acoustic model: ', acoustic_model ] );
end

