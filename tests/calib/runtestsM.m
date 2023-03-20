% calibration test cases
% mbp

% isovelocity calibration case

bellhopM( 'calibB' )
plotshd( 'calibB.shd.mat', 2, 2, 1 )
caxisrev( [ 40 80 ] );

bellhopM( 'calibB_gb' )
plotshd( 'calibB_gb.shd.mat', 2, 2, 2 )
caxisrev( [ 40 80 ] );

kraken( 'calibK' )
plotshd( 'calibK.shd.mat', 2, 2, 3 )
caxisrev( [ 40 80 ] );

scooterM( 'calibS' )
plotshd( 'calibS.shd.mat', 2, 2, 4 )
caxisrev( [ 40 80 ] );

%%
% gradient calibration case

bellhopM( 'calibBgrad' )
plotshd( 'calibBgrad.shd.mat', 2, 2, 1 )
caxisrev( [ 40 80 ] );

bellhopM( 'calibBgrad_gb' )
plotshd( 'calibBgrad_gb.shd.mat', 2, 2, 2 )
caxisrev( [ 40 80 ] );

kraken( 'calibKgrad' )
plotshd( 'calibKgrad.shd.mat', 2, 2, 3 )
caxisrev( [ 40 80 ] );

scooterM( 'calibSgrad' )
plotshd( 'calibSgrad.shd.mat', 2, 2, 4 )
caxisrev( [ 40 80 ] );
