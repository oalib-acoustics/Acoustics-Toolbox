clear nx ny Bathy

btyfil = 'provinces.bty';
interp_type = 'RL';

% Write a bathymetry from the workspace variables

% write bty file for BELLHOP3D
% mbp 3/2011

xctr = 30;
yctr = 30;

xmin = -10;
xmax = +10;
nx   = 11;

ymin = -10;
ymax = +10;
ny   = 21;

x = linspace( xmin, xmax, nx );
y = linspace( ymin, ymax, ny );

z = zeros( ny, nx );

for ix = 1 : nx
   for iy = 1 : ny
      z( iy, ix ) = 100.0;
   end
end

Bathy.X = x;
Bathy.Y = y;
Bathy.depth = z;

% set up matrix of bottom-type provinces
for ix = 1 : nx
   for iy = 1 : ny
      xt = x( ix );
      yt = y( iy );
      Bathy.province( iy, ix ) = 5;
      if ( xt >= 0.0 && yt >= 0 ) Bathy.province( iy, ix ) = 1; end
      if ( xt <  0.0 && yt >= 0 ) Bathy.province( iy, ix ) = 2; end
      if ( xt <  0.0 && yt <  0 ) Bathy.province( iy, ix ) = 3; end
      if ( xt >= 0.0 && yt <  0 ) Bathy.province( iy, ix ) = 4; end
      % if ( xt < 0.0 && yt > +xt + 5 ) Bathy.province( iy, ix ) = 1; end
      % if ( xt < 0.0 && yt < -xt - 5 ) Bathy.province( iy, ix ) = 2; end
      % if ( xt > 0.0 && yt > -xt + 5 ) Bathy.province( iy, ix ) = 3; end
      % if ( xt > 0.0 && yt < +xt - 5 ) Bathy.province( iy, ix ) = 4; end
   end
end

% set up definitions of bottom types

NProvinces = max( max( Bathy.province ) );

Bathy.geotype( 1, : ) = [ 1500 0 1.8 0.8 0 ];
Bathy.geotype( 2, : ) = [ 1550 0 1.8 0.8 0 ];
Bathy.geotype( 3, : ) = [ 1600 0 1.8 0.8 0 ];
Bathy.geotype( 4, : ) = [ 1650 0 1.8 0.8 0 ];
Bathy.geotype( 5, : ) = [ 1650 0 1.8 0.8 0 ];

writebdry3d( btyfil, interp_type, Bathy )


