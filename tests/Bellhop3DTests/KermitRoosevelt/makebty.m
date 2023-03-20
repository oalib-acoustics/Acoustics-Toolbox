btyfil = 'KR.bty';
interp_type = 'R';

nx = 251;
ny = 201;

% set up a flat bottom defined on a rectangular grid

Bathy.X = linspace( 0, 11000, nx );
Bathy.Y = linspace( -5000, 5000, ny );
Bathy.depth = 250 * ones( ny, nx );   % water depth of 250 m

% find x-y coordinates of the base of the seamount

cone.R = 350;   % radius in m
cone.H = 200;   % height in m
cone.x0 = 5000; % center x
cone.y0 = 0;    % center y

for ix = 1 : nx
    R = sqrt( ( Bathy.X( ix ) - cone.x0 )^2 + ( Bathy.Y - cone.y0 ).^2 );
    h = cone.H * ( 1 - R / cone.R );   % linear slope to cone height
    
    % indices covering the base of the cone
    iy = find( ( Bathy.X( ix ) - cone.x0 )^2 + ( Bathy.Y - cone.y0 ).^2 < cone.R.^2 );
    Bathy.depth( iy, ix ) = 250 - h( iy );
end

Bathy.X = Bathy.X / 1000;
Bathy.Y = Bathy.Y / 1000;

writebdry3d( btyfil, 'R', Bathy )

plotbdry3d( btyfil )

