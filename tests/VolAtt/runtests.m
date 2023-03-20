% Run tests to verify the volume attenuation is working correctly

global units
units = 'm';

% BELLHOP point source cases
bellhop 'freeB'
plotshd( 'freeB.shd', 3, 1, 1 );
axis( [ 0 5000 0 5000 ] )
caxisrev( [ 50 100 ] )

bellhop 'free_gbtB'
plotshd( 'free_gbtB.shd', 3, 1, 2 );
axis( [ 0 5000 0 5000 ] )
caxisrev( [ 50 100 ] )

% SCOOTER

scooter 'freeS'
plotshd( 'freeS.shd.mat', 3, 1, 3 );
axis( [ 0 5000 0 5000 ] )
caxisrev( [ 50 100 ] )

%%
% incoherent TL runs

% BELLHOP point source cases
bellhop 'freeB_Inc'
plotshd( 'freeB_Inc.shd', 3, 1, 1 );
axis( [ 0 5000 0 5000 ] )
caxisrev( [ 50 100 ] )

bellhop 'free_gbtB_Inc'
plotshd( 'free_gbtB_Inc.shd', 3, 1, 2 );
axis( [ 0 5000 0 5000 ] )
caxisrev( [ 50 100 ] )

plotshd( 'freeS.shd.mat', 3, 1, 3 );
axis( [ 0 5000 0 5000 ] )
caxisrev( [ 50 100 ] )

%%
% Thorp volume attenuation

% TL at 10 km would normally be 80 dB (spherical spreading)
% The volume attenuation adds about 3 dB

% BELLHOP point source cases
bellhop 'free_ThorpB'
plotshd( 'free_ThorpB.shd', 2, 1, 1 );
axis( [ 0 10000 0 5000 ] )
caxisrev( [ 50 100 ] )


% BELLHOP point source cases
bellhop 'free_FGB'
plotshd( 'free_FGB.shd', 2, 1, 2 );
axis( [ 0 10000 0 5000 ] )
caxisrev( [ 50 100 ] )

figure; plottlr( 'free_ThorpB.shd', 2500 )
hold on; plottlr( 'free_FGB.shd', 2500 )
%%
% SCOOTER point source cases

% This is commented out because it takes about 5 hours to run

% scooter 'free_ThorpS'
% plotshd( 'free_Thorp.shd.mat', 2, 1, 2 );
% axis( [ 0 10000 0 5000 ] )
% caxisrev( [ 50 100 ] )
% 
% scooter 'free_FGS'
% plotshd( 'free_FGS.shd.mat', 2, 1, 2 );
% axis( [ 0 10000 0 5000 ] )
% caxisrev( [ 50 100 ] )
% 
% figure; plottlr( 'free_ThorpS.shd', 2500 )
% hold on; plottlr( 'free_FGS.shd', 2500 )