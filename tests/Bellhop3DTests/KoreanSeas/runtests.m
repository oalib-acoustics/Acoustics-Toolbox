% run the Korean Sea 3D test cases
% mbp, July 2011

% the polar and sideview plots should be done in separate BELLHOP runs for
% speed ...

global units
units = 'km';

%%
%makebty;             % make the bathymetry

% Nx2D transmission loss

fileroot = 'KoreanSea_Nx2D'
   
copyfile( 'KoreanSea.bty', [ fileroot '.bty' ] )
copyfile( 'KoreanSea.ssp', [ fileroot '.ssp' ] )

bellhop3d( fileroot )     % run BELLHOP3D on the munk3d test case

% polar plot of the TL
figure
plotshdpol( [ fileroot '.shd' ], [ 950 1050 ], [ 750 900 ], 400 )
caxisrev( [ 60 100 ] )
axis equal

% side view of the TL
figure
plotshd(    [ fileroot '.shd' ] )
caxisrev( [ 60 100 ] )

delete( [ fileroot '.bty' ] )
delete( [ fileroot '.ssp' ] )

%% 3D transmission loss

fileroot = 'KoreanSea_3D'
   
copyfile( 'KoreanSea.bty', [ fileroot '.bty' ] )
copyfile( 'KoreanSea.ssp', [ fileroot '.ssp' ] )

bellhop3d( fileroot )     % run BELLHOP3D on the munk3d test case

% polar plot of the TL
figure
plotshdpol( [ fileroot '.shd' ], [ 950 1050 ], [ 750 900 ], 400 )
caxisrev( [ 60 100 ] )
axis equal

% side view of the TL
figure
plotshd(    [ fileroot '.shd' ] )
caxisrev( [ 60 100 ] )

delete( [ fileroot '.bty' ] )
delete( [ fileroot '.ssp' ] )


%% Nx2D ray trace

fileroot = 'KoreanSea_Nx2D_ray'

copyfile( 'KoreanSea.bty', [ fileroot '.bty' ] )
copyfile( 'KoreanSea.ssp', [ fileroot '.ssp' ] )

bellhop3d( fileroot )

figure
plotbdry3d( [ fileroot '.bty' ] )
hold on

plotray3d( [ fileroot '.ray' ] )
view( [ 70.5, 68 ] )

delete( [ fileroot '.bty' ] )
delete( [ fileroot '.ssp' ] )

%% 3D ray trace

fileroot = 'KoreanSea_3D_ray'

copyfile( 'KoreanSea.bty', [ fileroot '.bty' ] )
copyfile( 'KoreanSea.ssp', [ fileroot '.ssp' ] )

bellhop3d( fileroot )

figure
plotbdry3d( [ fileroot '.bty' ] )
hold on

plotray3d( [ fileroot '.ray' ] )
view( [ 70.5, 68 ] )

delete( [ fileroot '.bty' ] )
delete( [ fileroot '.ssp' ] )

%% 3D eigenrays

% This takes about 4 hours to run on my laptop ...
fileroot = 'KoreanSea_3D_eigen'

copyfile( 'KoreanSea.bty', [ fileroot '.bty' ] )
copyfile( 'KoreanSea.ssp', [ fileroot '.ssp' ] )

bellhop3d( fileroot )

figure
plotbdry3d( [ fileroot '.bty' ] )
hold on

plotray3d( [ fileroot '.ray' ] )
view( [ 70.5, 68 ] )

delete( [ fileroot '.bty' ] )
delete( [ fileroot '.ssp' ] )