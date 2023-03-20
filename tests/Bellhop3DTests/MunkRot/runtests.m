% run the 3D Munk case
% The angular limits were originally chosen so as not to hit the free surface
% check these closely if you see discrepancies

% Later on, I added walls in the bathymetry and these generate reflections

% The rotated case (Munk3D) is using walls to represent the top/bottom boundaries
% Since the bottom is all the same material we can't have different BCs on the two walls

% Fine sampling is required to represent vertical walls well

% 5/2019: I shifted the Munk SSP so that it goes from x= 0 to 5000 m
% The original version of BELLHOP3D required the source to be at the origin
% so I had shifted the SSP to the left to get the same effect as a source at 1000 m

global units
units = 'km';

%%
bellhop3d Munk3D

figure
plotshdpol( 'Munk3D.shd', 1, 0, 3000 )
caxisrev( [ 50 100 ] )
axis( [ 0 5 0 100 ] )

%%
bellhop3d Munk3Dz

figure
plotshd( 'Munk3Dz.shd' )
caxisrev( [ 50 100 ] )


%%
% The following test cases were to look very closely at the wall
% reflections to verify the treatment of the boundary conditions
bellhop3d foo

figure
plotshdpol( 'foo.shd', 0, 0, 3000 )
caxisrev( [ 50 100 ] )
axis( [ -1 4 0 100 ] )

%%
bellhop3d fooz

figure
plotshd( 'fooz.shd' )
caxisrev( [ 50 100 ] )

%%
bellhop3d foozz

figure
plotshd( 'foozz.shd' )
caxisrev( [ 50 100 ] )
