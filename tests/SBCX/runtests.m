% SBCX test cases
% These tests exercise the arrivals option in Bellhop
% verifying that the results obtained by summing from the arrivals file
% match those when BELLHOP directly outputs a pressure field.
%
% mbp

bellhop( 'sbcx' )
plotshd( [ 'sbcx.shd' ], 3, 1, 1 )
caxisrev( [ 40 80 ] )

bellhop( 'sbcx_Arr_asc' )
makeshdarr( 'sbcx_Arr_asc', 200, 'ascii' )
plotshd( 'sbcx_Arr_asc.shd.mat', 3, 1, 2 )
caxisrev( [ 40 80 ] )

bellhop( 'sbcx_Arr_bin' )
makeshdarr( 'sbcx_Arr_bin', 200, 'binary' )
plotshd( 'sbcx_Arr_bin.shd.mat', 3, 1, 3 )
caxisrev( [ 40 80 ] )


