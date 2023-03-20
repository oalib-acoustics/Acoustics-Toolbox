% generate bottom bathymetry producing a wall along y=4000


% mbp Feb. 15, 2019

x = linspace( -1, 6,  701 );   % needs fine sampling to get near vertical walls
y = linspace( -1, 101, 203 );
nx = length( x );
ny = length( y );

btyfil = 'Munk3D.bty';
atifil = 'Munk3D.ati';
interp_type = 'R';

[ X, Y ] = meshgrid( x, y );
z = 5000 * ones( ny, nx );

ii = find( X > 5 );
z( ii ) = 0.0;   % build a wall

ii = find( X < 0 );
z( ii ) = 0.0;   % build a wall

Bdry.X     = x;
Bdry.Y     = y;
Bdry.depth = z;
%%
writebdry3d( btyfil, interp_type, Bdry )

%%
%Bdry.depth = -z';
%writebdry3d( atifil, interp_type, Bdry )
