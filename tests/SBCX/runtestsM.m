% SBCX test cases
% These tests exercise the arrivals option in Bellhop
% verifying that the results obtained by summing from the arrivals file
% match those when BELLHOP directly outputs a pressure field.
%
% mbp

bellhopM( 'sbcx' )
figure
plotshd( [ 'sbcx.shd.mat' ] ) % , 2, 1, 1 )
caxis( [ 40 80 ] )

% following commented out because AddArr.m routine in BELLHOP is not
% implemented

% bellhopM( 'sbcx_Arr_asc' )
% makeshdarr( 'sbcx_Arr_asc', 200, 'matlab' )
% plotshd( 'sbcx_Arr_asc.mat', 2, 1, 2 )
% caxis( [ 40 80 ] )
