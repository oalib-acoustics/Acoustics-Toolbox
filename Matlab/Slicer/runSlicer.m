
% This script provides an example of how to run the environmental slicer

% There are 3 possible run types
%   'Grid'  TL on a regular 3D grid of receiver points
%   'Endpoints' TL on a slice between two endpoints (source and receiver)
%   'Polar' TL on a radial slice from the source

Params.RunType        = 'GRID';

Params.BathyFile      = 'Stellwagen.xyz';	% NGDC bathymetry for Stellwagen Bank
Params.acoustic_model = 'BELLHOP';              % either 'BELLHOP' or 'RAM'

Params.src_lat   =  42.55;                      % source latitude  in decimal degrees
Params.src_lon   = -70.55;	                % source longitude in decimal degrees
Params.src_depth = 15.0;                        % source depth in meters
Params.src_freq  = 100.0;                       % source frequency in Hz
Params.src_bw    = 0.075;                       % proportional bandwidth, [max(f)-min(f)] / mean(f)

%Params.Plane.src_lat     =  42.630;  % lat/lon of the source
%Params.Plane.src_lon     = -70.600;

switch Params.RunType
   case 'GRID'
      % Set the lat, lon, depth grid where TL is interpolated.
      %
      % Each is specified by the start and end values and the number of grid points.
      % Southern latitudes and Western longitudes are negative (decimal degrees).
      
      Params.grid_lat = [ 42.5,  42.8, 251 ];
      Params.grid_lon = [-70.8, -70.2, 201 ];
      Params.grid_depths = [ 10, 50, 5 ];	% depths: 10, 20, 30, 40, 50 m
      
      % The following parameters control the resolution of the underlying polar
      % grid where TL is calculated.
      % This is a fan of vertical planes. Using lower resolution increases the
      % accuracy of the interpolation (at the expense of more computer time).
      
      Params.bearing_resolution  = 2.0;	% resolution of polar angles (degrees)
      Params.range_resolution_km = 0.05;	% resolution of range receivers (km)
      
   case 'ENDPOINTS'
      
      % Set the receiver geometry of the computational plane / slice
      
      Params.rcv_lat    =  42.720;  % lat/lon of most distant receivers
      Params.rcv_lon    = -70.599;
      
      Params.range_resolution_km = 0.05;	% resolution of range receivers (km)
      Params.grid_depths = [ 0, 300, 201 ];	% depths: 10, 20, 30, 40, 50 m
      
   case 'POLAR'
      
      % receiver geometry of the computational plane / slice
      
      Params.bearing_deg = 2.0;      % bearing of receiver plane (deg cw from true N)
      Params.max_r_km    = 10;       % range of most distant receivers in km
      
      Params.range_resolution_km = 0.05;	% resolution of range receivers (km)
      Params.grid_depths = [ 0, 300, 201 ];	% depths: 10, 20, 30, 40, 50 m
      
end

Params.SSP_month = 6;	% month of year for WOA SSP, 1=Jan, 2=Feb ...

% Sediments in Stellwagen Bank are generally Terrigenous Gravel, Sand.
%
% geo-acoustic parameters for grain size of phi = 1, coarse/med sand (Wentworth)
% Params.Bot.alphaR = 1806.0;
% Params.Bot.betaR  = 0.0;
% Params.Bot.rho    = 2.151;
% Params.Bot.alphaI = 0.497;
% Params.Bot.betaI  = 0.0;
%
% geo-acoustic parameters for grain size of phi = 2, med/fine sand (Wentworth)
% Params.Bot.alphaR = 1681.0;
% Params.Bot.betaR  = 0.0;
% Params.Bot.rho    = 1.615;
% Params.Bot.alphaI = 0.523;
% Params.Bot.betaI  = 0.0;
%
% geo-acoustic parameters for grain size of phi = 3, fine/vf sand (Wentworth)
Params.Bot.alphaR = 1593.0;
Params.Bot.betaR  = 0.0;
Params.Bot.rho    = 1.339;
Params.Bot.alphaI = 0.592;
Params.Bot.betaI  = 0.0;
%
% geo-acoustic parameters for grain size of phi = 4, vf sand/silt (Wentworth)
% Params.Bot.alphaR = 1529.0;
% Params.Bot.betaR  = 0.0;
% Params.Bot.rho    = 1.224;
% Params.Bot.alphaI = 0.721;
% Params.Bot.betaI  = 0.0;

% Run the tl_grid function, it returns TLmat( nlat, nlon, ndepths )

TLmat = Slicer( Params );
%%

% Example plot

figure
tej = flipud( jet( 256 ) );  % 'jet' colormap reversed
colormap( tej )

switch Params.RunType
   case 'GRID'
      lat_deg = linspace( Params.grid_lat( 1 ), Params.grid_lat( 2 ), Params.grid_lat( 3 ) );
      lon_deg = linspace( Params.grid_lon( 1 ), Params.grid_lon( 2 ), Params.grid_lon( 3 ) );
      
      imagesc( lon_deg, lat_deg, TLmat( :, :, 2 ) )
      axis xy
      title( 'Stellwagen Bank Example - Receiver depth = 20 m' );
      xlabel( 'Longitude (degrees)' );
      ylabel( 'Latitude (degrees)' );
      
      caxisrev( [ 50 80 ] )
      
   case 'POLAR'
      depth = linspace( Params.grid_depths( 1 ), Params.grid_depths( 2 ), Params.grid_depths( 3 ) );
      range = linspace( 0, Params.max_r_km, 50 );
      
      imagesc( range, depth, TLmat )
      colorbar
      xlabel( 'Range (km)' )
      ylabel( 'Depth (m)' )
      
      caxisrev( [ 50 80 ] )
     
   case 'ENDPOINTS'
      lat_deg = linspace( Params.src_lat, Params.rcv_lat, 50 );
      lon_deg = linspace( Params.src_lon, Params.rcv_lon, 50 );
      
      imagesc( lon_deg, depth, TLmat )
 
      colorbar
      xlabel( 'Longitude (degrees)' )
      ylabel( 'Depth (m)' )
      
      caxisrev( [ 50 80 ] )
end
