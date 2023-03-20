function TLmat = Slicer( Params )

% TL_GRID - compute transmission loss on a rectangular grid (in lat/lon)

persistent BathyData
persistent SlicerParams
persistent SSPmonth SSPdata

% List of the WOA monthly files

WOA_files = { '01-January.mat', ...
   '02-February.mat', ...
   '03-March.mat', ...
   '04-April.mat', ...
   '05-May.mat', ...
   '06-June.mat', ...
   '07-July.mat', ...
   '08-August.mat', ...
   '09-September.mat', ...
   '10-October.mat', ...
   '11-November.mat', ...
   '12-December.mat' };

Slice = CreateSliceStructure( );   % Create a new Slice struct with default values

% following should only be done on the first call
[ SlicerParams, BathyData ] = slicer_init;

% Get a sound speed profile from the WOA

if isempty( SSPmonth )
   % persistent variable has not been initialized
   SSPmonth = Params.SSP_month;
end

if isempty( SSPdata ) || ~strcmp( SSPmonth, Params.SSP_month )
   % input value has changed, reload the sound speed profile data
   SSPdata = load( [ 'WOA/', WOA_files{SSPmonth} ] );
end

%%

src_lat = Params.src_lat;
src_lon = Params.src_lon;

Ndepths = length( SSPdata.Oceanography.depth );
[ ~, i_lat ] = min( abs( SSPdata.Oceanography.lat - src_lat ) );
[ ~, i_lon ] = min( abs( SSPdata.Oceanography.lon - src_lon ) );

% Populate the Slice struct with the user's values

Slice.acoustic_model = Params.acoustic_model;

Slice.freq = Params.src_freq;

Slice.Beam.RunType( 1 ) = 'C';

Slice.Bdry.Bot.alphaR = Params.Bot.alphaR;
Slice.Bdry.Bot.betaR  = Params.Bot.betaR;
Slice.Bdry.Bot.rho    = Params.Bot.rho;
Slice.Bdry.Bot.alphaI = Params.Bot.alphaI;
Slice.Bdry.Bot.betaI  = Params.Bot.betaI;

Slice.Pos.s.z = Params.src_depth;

Slice.SSP.depth(2)   = SSPdata.Oceanography.depth( end );
Slice.SSP.raw.z      = SSPdata.Oceanography.depth( : );
Slice.SSP.raw.alphaR = squeeze( SSPdata.Oceanography.c( i_lon, i_lat, : ) );
Slice.SSP.raw.betaR  = zeros( 1, Ndepths );
Slice.SSP.raw.rho    =  ones( 1, Ndepths );
Slice.SSP.raw.alphaI = zeros( 1, Ndepths );
Slice.SSP.raw.betaI  = zeros( 1, Ndepths );

switch Params.RunType
   case 'GRID'
      %%
      % Setup the lat, lon grid where TL values will be computed
      
      grid_lat = Params.grid_lat;
      grid_lon = Params.grid_lon;
      
      nlat = grid_lat( 3 );
      nlon = grid_lon( 3 );
      
      lat = linspace( grid_lat( 1 ), grid_lat( 2 ), nlat );
      lon = linspace( grid_lon( 1 ), grid_lon( 2 ), nlon );
      %%
      
      % Find a suitable branch cut for the bearing angles of the rectangular box
      % for the receiver grid
      branch = [ 0 360 ];   % default branch is a full sweep
      
      % get bearings to the 4 corners of the receiver box
      [ ~, bearing1, ~ ] = dist_wh_math( [ src_lat grid_lat( 1 ) ], [ src_lon grid_lon( 1 ) ] );
      [ ~, bearing2, ~ ] = dist_wh_math( [ src_lat grid_lat( 1 ) ], [ src_lon grid_lon( 2 ) ] );
      [ ~, bearing3, ~ ] = dist_wh_math( [ src_lat grid_lat( 2 ) ], [ src_lon grid_lon( 2 ) ] );
      [ ~, bearing4, ~ ] = dist_wh_math( [ src_lat grid_lat( 2 ) ], [ src_lon grid_lon( 1 ) ] );
      
      if     bearing1 >=  90 && bearing1 < 180   % Q2
         branch( 2 ) = bearing1;
      elseif bearing1 >= 270 && bearing1 < 360   % Q4
         branch( 1 ) = bearing1;
      end
      
      if     bearing2 >=   0 && bearing2 <  90   % Q1
         branch( 1 ) = bearing2;
      elseif bearing2 >= 180 && bearing2 < 270   % Q3
         branch( 2 ) = bearing2;
      end
      
      if     bearing3 >=  90 && bearing3 < 180   % Q2
         branch( 1 ) = bearing3;
      elseif bearing3 >= 270 && bearing3 < 360   % Q4
         branch( 2 ) = bearing3;
      end
      
      if     bearing4 >=   0 && bearing4 <  90   % Q1
         branch( 2 ) = bearing4;
      elseif bearing4 >= 180 && bearing4 < 270   % Q3
         branch( 1 ) = bearing4;
      end
      
      % Compute the bearing and range values for all lat, lon grid points
      
      grid_range_km    = zeros( nlat, nlon );
      grid_bearing_deg = zeros( nlat, nlon );
      range_m          = zeros( nlat, 1    );
      bearing_12       = zeros( nlat, 1    );
      
      % following should be vectorized (MPB)
      for j = 1 : nlon
         
         for i = 1 : nlat
            % compute range and bearing from source to this grid point
            [ range_m( i ), bearing_12( i ), ~ ] = dist_wh_math( [ src_lat, lat( i ) ], ...
               [ src_lon, lon( j ) ] );
         end
         
         % adjust the bearing to be consistent with the branch cut used
         ii = find( bearing_12 < branch( 1 ) );
         bearing_12( ii ) = bearing_12( ii ) + 360;
         
         ii = find( bearing_12 > branch( 2 ) );
         bearing_12( ii ) = bearing_12( ii ) - 360;
         
         % save the bearing and range
         grid_range_km( :, j )    = range_m / 1000.0;
         grid_bearing_deg( :, j ) = bearing_12;
      end
      %%
      
      % Determine the min, max values of both the bearing angles and range
      
      bearing_min = min( min( grid_bearing_deg ) );
      bearing_max = max( max( grid_bearing_deg ) );
      
      % Apply any needed tweaks to the min, max of the bearing angles
      
      if ( bearing_max - bearing_min ) > 180
         % calculate over a full disk
         bearing_min = branch( 1 );
         bearing_max = branch( 2 );
      end
      
      nbearings = round( ( bearing_max - bearing_min ) / Params.bearing_resolution );
      nbearings = max( [ nbearings, 2 ] );
      
      bearings = linspace( bearing_min, bearing_max, nbearings );
      
      % Apply any needed tweaks to the min, max of the range
      
      range_min   = min( min( grid_range_km ) );
      range_max   = max( max( grid_range_km ) );
      
   case 'ENDPOINTS'
      nbearings = 1;
      
      rcv_lat = Params.rcv_lat;
      rcv_lon = Params.rcv_lon;
      
      % Compute the bearing and range of the receiver (relative to src)
      
      [ r_m, bearing_deg, ~ ] = dist_wh_math( [ src_lat rcv_lat ], [ src_lon rcv_lon ] );
      bearings( 1 ) = bearing_deg;      % bearing of receiver plane (deg cw from true N)

      range_min   = 0;
      range_max   = Params.max_r_km;
      
   case 'POLAR'
      
      nbearings = 1;
      bearings( 1 ) = Params.bearing_deg;      % bearing of receiver plane (deg cw from true N)
      
      range_min   = 0;
      range_max   = Params.max_r_km;
      
   otherwise
      error( [ mfilename, ': unrecognized method for specification of computational plane: ', method ] );
end


range_resolution_km = Params.range_resolution_km;

if Params.src_bw > 0.0
   alpha = Params.src_bw;
   
   % case of range averaging (as a surrogate for frequency averaging)
   range_min = ( 1.0 - alpha ) * range_min;
   range_max = ( 1.0 + alpha ) * range_max;
   del_range = range_max - range_min;
   
   % range averaging also influences the (required) range resolution
   ra_resolution_km  = 0.5 * alpha * ( 0.9 * range_min + 0.1 * range_max ) / 3.0;
   min_resolution_km = min( [ range_resolution_km, ra_resolution_km ] );
   nrr               = 1 + round( del_range / min_resolution_km );
else
   % no range averaging computation
   range_min = max( [ range_min - 0.20, 0.0 ] );
   range_max = range_max + 0.20;
   del_range = range_max - range_min;
   nrr       = 1 + round( del_range / range_resolution_km );
end

% nrr       = max( [ nrr, nlat, nlon, 2 ] );
nrr       = max( [ nrr, 2 ] );

%%

% Finish setup of the Slice struct

grid_depths = Params.grid_depths;

min_d_m = grid_depths( 1 );
max_d_m = grid_depths( 2 );
nrd     = grid_depths( 3 );

Plane.method      = 'POLAR';
Plane.src_lat     = src_lat;
Plane.src_lon     = src_lon;
Plane.min_r_km    = range_min;
Plane.max_r_km    = range_max;
Plane.min_d_m     = min_d_m;
Plane.max_d_m     = max_d_m;
Plane.nrr         = nrr;
Plane.nrd         = nrd;

Slice.Plane = Plane;

%%

% Loop over bearings and calculate TL

TLpolar = zeros( nrd, nrr, nbearings );

for iangle = 1 : nbearings
   
   disp( [ 'iangle = ' num2str( iangle ), ' of ', num2str( nbearings ) ] )
   
   Slice.Plane.bearing_deg = bearings( iangle );
   
   % *** Run the acoustic propagation model on the Slice ***
   
   [ Slice_out, tl_db, rc ] = SingleSlice( SlicerParams, Slice, BathyData );
   
   if rc ~= 0   % run went OK?
      error( [ mfilename, ': execution of acoustic model (', Slice.acoustic_model, ...
         ') returned error code: ', num2str( rc ) ] );
   end
   
   % if first bearing, pre-compute window for range averaging
   if iangle == 1
      % save the receiver ranges
      r_km  = Slice_out.Pos.r.r;
      dr_km = r_km( 2 ) - r_km( 1 );
      
      % One time computations when range averaging is enabled
      if Params.src_bw > 0.0
         % number of points in the half-window for range averaging
         no_pts = round( ( 0.5 * alpha * r_km ) ./ dr_km );
         no_pts = max( no_pts, 1 );   % make sure a minimum of 1 point in the average
         
         % indices where range averaging computation is feasible
         iranges = find( ( ( 1 : nrr ) - no_pts + 1 ) >= 1 & ...
            ( ( 1 : nrr ) + no_pts - 1 ) <= nrr );
      end
   end
   
   % Perform range averaging to approximate averaging over acoustic frequency
   
   if Params.src_bw > 0.0
      tl = 10.0 .^ ( tl_db / -20.0 );  % convert TL from dB to linear units
      
      % compute range averages for each depth
      tl_bar = tl;
      
      for ir = iranges
         lim1 = ir - no_pts( ir ) + 1;
         lim2 = ir + no_pts( ir ) - 1;
         tl_bar( :, ir ) = mean( tl( :, ( lim1 : lim2 ) ), 2 );
      end
      
      tl_db = -20.0 * log10( tl_bar );   % convert the range averaged TL to dB
   end
   
   % need to specify 1 : nrr below because RAM may return an extra point
   TLpolar( :, :, iangle ) = tl_db( :, 1 : nrr );   % Store this slice for later spatial interpolation
   
end

%%
% Return TL result in appropriate form

switch Params.RunType
   case 'GRID'
      % For each depth, interpolate from the polar grid to the caller's lat, lon grid
      
      TLmat = zeros( nlat, nlon, nrd );
      
      for idepth = 1 : nrd
         % interpolate to the exact bearing and range of this lat, lon point
         TLmat( :, :, idepth ) = interp2( bearings, r_km, ...
            squeeze( TLpolar( idepth, :, : ) ), ...
            grid_bearing_deg, grid_range_km );
      end
      
   case { 'ENDPOINTS', 'POLAR' }
      TLmat = TLpolar( :, :, : );
end
