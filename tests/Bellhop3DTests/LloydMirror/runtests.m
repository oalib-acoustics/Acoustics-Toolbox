
% Test of BELLHOP3D surface reflection (Lloyd mirror). This test is
% just the 3D version of the corresponding test for vanilla BELLHOP.

bellhop3d 'LloydSurfaceB3'
plotshd( 'LloydSurfaceB3.shd', 2, 2, 1 );
caxisrev( [ 20 90 ] )
colorbar( 'YDir', 'Reverse' )

% exact solution
freq = 5;
k0 = 2 * pi * freq / 1500.0;
zs = 2500;
zr =  0 : 25 : 5000;
rr = 25 : 25 : 10000;

ppoint = zeros( length( zr ), length( rr )  );

for iz = 1: length( zr )
   rmat1 = sqrt( rr.^2 + ( zr( iz ) - zs )^2 );
   rmat2 = sqrt( rr.^2 + ( zr( iz ) + zs )^2 );
   ppoint( iz, : ) = exp( 1i * k0 * rmat1 ) ./ rmat1 - exp( 1i * k0 * rmat2 ) ./ rmat2;
end

subplot( 2, 2, 2 )
pcolor( rr, zr, -20 * log10( abs( ppoint ) ) );
shading flat;
colorbar( )
set( gca, 'YDir', 'Reverse' )
set( gca, 'TickDir', 'out' )
set( findall( gcf, 'type', 'ColorBar' ), 'TickDir', 'out' )
caxisrev( [ 20 90 ] )
title( sprintf('Analytical solution, surface reflection\nFreq = %d Hz    z_{src} = %d m', freq, zs) )
xlabel( 'Range (m)' )
ylabel( 'Depth (m)' )

% BELLHOP3D wall reflection (Lloyd mirror) test. This test is similar
% to the one above, but for reflection from a nearly vertical wall. It
% compares BELLHOP3D in full 3D mode to vanilla BELLHOP. There are some
% small differences in the depth of the nulls. This is to be expected
% as the BELLHOP3D wall is flat, while the BELLHOP wall is actually
% cylindrical which produces some additional focusing of the reflected
% rays. The intent is to verify the phasing of the reflected wave.

bellhop3d 'LloydWall_1B3'
plotshd( 'LloydWall_1B3.shd', 2, 2, 3 );
caxisrev( [ 20 90 ] )
colorbar( 'YDir', 'Reverse' )

bellhop 'LloydWallB'
plotshd( 'LloydWallB.shd', 2, 2, 4 );
caxisrev( [ 20 90 ] )
colorbar( 'YDir', 'Reverse' )

% BELLHOP3D wall reflection (Lloyd mirror) test. This is similar to
% the one above. It compares walls at x = +/- 5 km, y = +/- 5 km. They
% should be the same.

bellhop3d 'LloydWall_1B3'
plotshd( 'LloydWall_1B3.shd', 2, 2, 1 );
caxisrev( [ 20 90 ] )
colorbar( 'YDir', 'Reverse' )

bellhop3d 'LloydWall_2B3'
plotshd( 'LloydWall_2B3.shd', 2, 2, 2 );
caxisrev( [ 20 90 ] )
colorbar( 'YDir', 'Reverse' )

bellhop3d 'LloydWall_3B3'
plotshd( 'LloydWall_3B3.shd', 2, 2, 3 );
caxisrev( [ 20 90 ] )
colorbar( 'YDir', 'Reverse' )

bellhop3d 'LloydWall_4B3'
plotshd( 'LloydWall_4B3.shd', 2, 2, 4 );
caxisrev( [ 20 90 ] )
colorbar( 'YDir', 'Reverse' )

%
% End of runtests.m
