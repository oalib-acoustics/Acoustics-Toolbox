% Run tests to verify the BELLHOP option for an irregular (vs. rectangular)
% set of receiver coordinates.
% The main use of this is to have a terrain following bottom receiver,
% which in turn is needed for reverb caculations

units = 'm';

bellhop( 'lower_half' )
plotshd( 'lower_half.shd', 2, 1, 1 );
axis( [ 0 10000 60 110 ] )

bellhop( 'lower_half_gbt' )
plotshd( 'lower_half_gbt.shd', 2, 1, 2 );
axis( [ 0 10000 60 110 ] )