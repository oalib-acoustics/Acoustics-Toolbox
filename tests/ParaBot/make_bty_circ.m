% generate parabolic bottom bathymetry
% mbp Dec. 27, 2005
% example from:
% R. W. McGirr, D.B. King, J.A. Davis, J. Campbell, "An evaluation of
% range-dependent ray theory models", NORDA report 115.

b = 250000;
c = 250;

r = 0:100:24750;
z = 0.002 * b * sqrt( 1 + r / c );

clear r z
z = 500:10:5000;
r = sqrt( 5000^2 - ( z - 5000 ).^2 );

fid = fopen( 'ParaBot.bty', 'w' );
fprintf( fid, '''%c'' \n', 'C' );
fprintf( fid, '%i', length( r ) );
for ii = 1: length( r )
   fprintf( fid, '\n %f %f', r( ii )/1000, z( ii ) );
end
fclose( fid );
