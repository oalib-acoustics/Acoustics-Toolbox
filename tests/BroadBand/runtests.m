% Testing the power law attenuation in sediments and the broadband option

% The Munk case verifies the broadband option with both KRAKEN and SCOOTER
% We use a big frequency spread largely as a check that the SCOOTER
% stabilising attenuation has been scaled properly with frequency (in both
% SCOOTER and FIELDSCO)

% We do a separate run at a single frequency (not broadband) to verify that
% the broadband version replicates the original narrowband one
% KRAKEN and SCOOTER results are virtually identical

% mbp
global units
units = 'km';

%%
kraken(  'MunkK' )

if ( exist( 'MunkK.shd.mat', 'file' ) )
   plotshd( 'MunkK.shd.mat',   50, 2, 1, 1 )
   plotshd( 'MunkK.shd.mat',  500, 2, 1, 2 )   
   %caxisrev( [ 50 100 ] )
end

%%
scooter( 'MunkS' )

if ( exist( 'MunkS.shd.mat', 'file' ) )
   plotshd( 'MunkS.shd.mat',   50, 3, 1, 1 )
   plotshd( 'MunkS.shd.mat',  500, 3, 1, 2 )
   %caxisrev( [ 50 100 ] )
end

%%
scooter( 'MunkS_500Hz' )

if ( exist( 'MunkS_500Hz.shd.mat', 'file' ) )
   plotshd( 'MunkS_500Hz.shd.mat', 500, 3, 1, 3 )
   %caxisrev( [ 50 100 ] )
end
