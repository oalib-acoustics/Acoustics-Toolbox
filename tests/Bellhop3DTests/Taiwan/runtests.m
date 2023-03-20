% run the Taiwan 3D test cases
% mbp, July 2011

% the polar and sideview plots should be done in separate BELLHOP runs for
% speed ...

global units
units = 'km';

%Bathy = LoadBathymetry( 'Taiwan.xyz' );   % load the xyz file from ETOPO1
%writebdry3d( 'Taiwan.bty', Bathy );        % save it in my bty format
figure
plotbdry3d Taiwan.bty                          % howzit look?
caxis( [ -7500 7500 ] )

%parula_dat = parula;
jet_dat = jet( 60 );

%ocean = flipud( parula_dat( 1 : 30, : ));
ocean = flipud( jet_dat( 1 : 30, : ));

land_dat  = copper( 60 );

ocean_land = [ land_dat( 15 : 44, : ) ; ocean ];
colormap( ocean_land )

%%
%figure
%plotssp3d Taiwan.ssp

% makebty;             % make the bathymetry
%%

% 2D Transmission Loss runs

fileroot = 'TaiwanNx2D'
   
copyfile( 'Taiwan.bty', [ fileroot '.bty' ] )
copyfile( 'Taiwan.ssp', [ fileroot '.ssp' ] )

bellhop3d( fileroot )

% polar plot of the TL
figure
plotshdpol( [ fileroot '.shd' ], [ 175 300 ], [ 200 300 400 ], 500 )
caxisrev( [ 60 100 ] )

% sideview of the TL
figure
plotshd(    [ fileroot '.shd' ] )
caxisrev( [ 60 100 ] )

delete( [ fileroot '.bty' ] )
delete( [ fileroot '.ssp' ] )

%%

% 3D Transmission Loss runs

fileroot = 'Taiwan3D'
   
copyfile( 'Taiwan.bty', [ fileroot '.bty' ] )
copyfile( 'Taiwan.ssp', [ fileroot '.ssp' ] )

bellhop3d( fileroot )

% polar plot of the TL
figure
plotshdpol( [ fileroot '.shd' ], [ 175 300 ], [ 200 300 400 ], 500 )
caxisrev( [ 60 100 ] )

print -depsc2 Taiwan3DTL

% sideview of the TL
figure
plotshd(    [ fileroot '.shd' ] )
caxisrev( [ 60 100 ] )

print -depsc2 Taiwan3D_side

delete( [ fileroot '.bty' ] )
delete( [ fileroot '.ssp' ] )


%%

% ray trace runs

fileroot = 'TaiwanNx2D_ray'

copyfile( 'Taiwan.bty', [ fileroot '.bty' ] )
copyfile( 'Taiwan.ssp', [ fileroot '.ssp' ] )

bellhop3d( fileroot )

figure
plotbdry3d( [ fileroot '.bty' ] )
caxis( [ -7500 7500 ] )
colormap( ocean_land )

hold on

plotray3d( [ fileroot '.ray' ] )

print -depsc2 TaiwanNx2D_ray

delete( [ fileroot '.bty' ] )
delete( [ fileroot '.ssp' ] )

%%

fileroot = 'Taiwan3D_ray'

copyfile( 'Taiwan.bty', [ fileroot '.bty' ] )
copyfile( 'Taiwan.ssp', [ fileroot '.ssp' ] )

bellhop3d( fileroot )

figure
plotbdry3d( [ fileroot '.bty' ] )
caxis( [ -7500 7500 ] )
colormap( ocean_land )

hold on

plotray3d( [ fileroot '.ray' ] )
axis( [ 0 350 150 600 -7500 7500 ] )

print -depsc2 Taiwan3D_ray

%%
delete( [ fileroot '.bty' ] )
delete( [ fileroot '.ssp' ] )
