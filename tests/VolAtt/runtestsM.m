% Run tests to verify the volume attenuation is working correctly

global units
units = 'm';

% BELLHOP point source cases
bellhopM 'freeB'
plotshd( 'freeB.shd.mat', 3, 1, 1 );
axis( [ 0 5000 0 5000 ] )
caxisrev( [ 50 100 ] )

bellhopM 'free_gbtB'
plotshd( 'free_gbtB.shd.mat', 3, 1, 2 );
axis( [ 0 5000 0 5000 ] )
caxisrev( [ 50 100 ] )

% SCOOTER

scooterM 'freeS'
plotshd( 'freeS.shd.mat', 3, 1, 3 );
axis( [ 0 5000 0 5000 ] )
caxisrev( [ 50 100 ] )

%%
% incoherent TL runs

% BELLHOP point source cases
bellhopM 'freeB_Inc'
plotshd( 'freeB_Inc.shd.mat', 3, 1, 1 );
axis( [ 0 5000 0 5000 ] )
caxisrev( [ 50 100 ] )

bellhopM 'free_gbtB_Inc'
plotshd( 'free_gbtB_Inc.shd.mat', 3, 1, 2 );
axis( [ 0 5000 0 5000 ] )
caxisrev( [ 50 100 ] )

plotshd( 'freeS.shd.mat', 3, 1, 3 );
axis( [ 0 5000 0 5000 ] )
caxisrev( [ 50 100 ] )

%%
% Thorp volume attenuation

% BELLHOP point source cases
bellhopM 'free_ThorpB'
plotshd( 'free_ThorpB.shd.mat', 2, 1, 1 );
axis( [ 0 10000 0 5000 ] )
caxisrev( [ 50 100 ] )


% BELLHOP point source cases
bellhopM 'free_FGB'
plotshd( 'free_FGB.shd.mat', 2, 1, 2 );
axis( [ 0 10000 0 5000 ] )
caxisrev( [ 50 100 ] )

figure; plottlr( 'free_ThorpB.shd.mat', 2500 )
hold on; plottlr( 'free_FGB.shd.mat', 2500 )