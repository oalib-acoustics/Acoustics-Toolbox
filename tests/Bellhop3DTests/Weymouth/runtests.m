% run the Weymouth 3D test cases
% mbp, July 2011

% the polar and sideview plots should be done in separate BELLHOP runs for
% speed ...

global units
units = 'm';

%Bathy = LoadBathymetry( 'Weymouth.xyz' );   % load the xyz file from ETOPO1
%writebdry3d( 'Weymouth.bty', Bathy );        % save it in my bty format
figure
plotbdry3d Weymouth.bty                          % howzit look?

% makebty;             % make the bathymetry
%%

% ray trace runs

fileroot = 'WeymouthNx2D_ray'

copyfile( 'Weymouth.bty', [ fileroot '.bty' ] )

bellhop3d( fileroot )

figure
plotbdry3d( [ fileroot '.bty' ] )
hold on

plotray3d( [ fileroot '.ray' ] )

delete( [ fileroot '.bty' ] )

print -dpng ray_Nx2D

%%


fileroot = 'Weymouth3D_ray'

copyfile( 'Weymouth.bty', [ fileroot '.bty' ] )

bellhop3d( fileroot )

figure
plotbdry3d( [ fileroot '.bty' ] )
hold on

plotray3d( [ fileroot '.ray' ] )

delete( [ fileroot '.bty' ] )
print -dpng ray_3D
%%

% 2D Transmission Loss runs

fileroot = 'WeymouthNx2D'
   
copyfile( 'Weymouth.bty', [ fileroot '.bty' ] )

bellhop3d( fileroot )

% polar plot of the TL
figure
%plotshdpol( [ fileroot '.shd' ], [ 0 ], [ -.25 ], 10 )
plotshdpol( [ fileroot '.shd' ], [ -.1 ], [ -.25 ], 10 )
caxisrev( [ 10 50 ] )
%axis( [ -.35 .35 -.4 .25 ] )
print -dpng TL_Nx2D

% sideview of the TL
figure
plotshd(    [ fileroot '.shd' ] )
caxisrev( [ 10 50 ] )

delete( [ fileroot '.bty' ] )

print -dpng TLside_Nx2D

%%

% 3D Transmission Loss runs

fileroot = 'Weymouth3D'
   
copyfile( 'Weymouth.bty', [ fileroot '.bty' ] )

bellhop3d( fileroot )

% polar plot of the TL
figure
%plotshdpol( [ fileroot '.shd' ], [ 0 ], [ -.25 ], 10 )
plotshdpol( [ fileroot '.shd' ], [ -.1 ], [ -.25 ], 10 )
caxisrev( [ 10 50 ] )
%axis( [ -.35 .35 -.4 .25 ] )
print -dpng TL_3D

% sideview of the TL
figure
plotshd(    [ fileroot '.shd' ] )
caxisrev( [ 10 50 ] )

delete( [ fileroot '.bty' ] )

print -dpng TLside_3D
