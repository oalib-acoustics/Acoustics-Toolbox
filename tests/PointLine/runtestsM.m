% Run tests to verify the point vs. line source options

% Geometric hat beams

% point source

bellhopM 'LloydPointB'
plotshd( 'LloydPointB.shd.mat', 2, 2, 1 );
caxisrev( [ 40 80 ] )
colorbar( 'YDir', 'Reverse' )

% line source

bellhopM 'LloydLineB'
plotshd( 'LloydLineB.shd.mat', 2, 2, 2 );
caxisrev( [ 0 25 ] )
colorbar( 'YDir', 'Reverse' )


%%

% exact solution
freq = 5;
k0 = 2 * pi * freq / 1500.0;
zs = 3000;
zr = 0:25:5000;
rr = 25:25:10000;
for iz = 1: length( zr )
   rmat1 = sqrt( rr.^2 + (zr( iz ) - zs )^2 );
   rmat2 = sqrt( rr.^2 + (zr( iz ) + zs )^2 );
   ppoint( iz, : ) = exp( i * k0 * rmat1 ) ./ rmat1 - exp( i * k0 * rmat2 ) ./ rmat2;
   pline(  iz, : ) = besselh( 0, k0 * rmat1 ) - besselh( 0, k0 * rmat2 );
end

pline = sqrt( pi / 2 ) * pline; %normalization to roughly match acoustics toolbox

subplot( 2, 2, 3 )
pcolor( rr, zr, -20 * log10( abs( ppoint ) ) ); shading interp; view( 0, -90 )
caxisrev( [ 40 80 ] )
colorbar( 'YDir', 'Reverse' )
title( 'Exact point source solution' )

subplot( 2, 2, 4 )
pcolor( rr, zr, -20 * log10( abs( pline ) ) ); shading interp; view( 0, -90 )
caxisrev( [  0 25 ] )
colorbar( 'YDir', 'Reverse' )
title( 'Exact line source solution' )

%%

% Geometric Gaussian beams

% point source

bellhopM 'LloydPoint_gbtB'
plotshd( 'LloydPoint_gbtB.shd.mat', 2, 2, 1 );
caxisrev( [ 40 80 ] )
colorbar( 'YDir', 'Reverse' )

% line source

bellhopM 'LloydLine_gbtB'
plotshd( 'LloydLine_gbtB.shd.mat', 2, 2, 2 );
caxisrev( [ 0 25 ] )
colorbar( 'YDir', 'Reverse' )

subplot( 2, 2, 3 )
pcolor( rr, zr, -20 * log10( abs( ppoint ) ) ); shading interp; view( 0, -90 )
caxisrev( [ 40 80 ] )
colorbar( 'YDir', 'Reverse' )
title( 'Exact point source solution' )

subplot( 2, 2, 4 )
pcolor( rr, zr, -20 * log10( abs( pline ) ) ); shading interp; view( 0, -90 )
caxisrev( [  0 25 ] )
colorbar( 'YDir', 'Reverse' )
title( 'Exact line source solution' )


%%

scooterM 'LloydSPoint'
plotshd( 'LloydSPoint.shd.mat', 2, 2, 1 );
caxisrev( [ 40 80 ] )
colorbar( 'YDir', 'Reverse' )

scooterM 'LloydSLine'
plotshd( 'LloydSLine.shd.mat', 2, 2, 2 );
caxisrev( [ 0 25 ] )
colorbar( 'YDir', 'Reverse' )

subplot( 2, 2, 3 )
pcolor( rr, zr, -20 * log10( abs( ppoint ) ) ); shading interp; view( 0, -90 )
caxisrev( [ 40 80 ] )
colorbar( 'YDir', 'Reverse' )
title( 'Exact point source solution' )

subplot( 2, 2, 4 )
pcolor( rr, zr, -20 * log10( abs( pline ) ) ); shading interp; view( 0, -90 )
caxisrev( [  0 25 ] )
colorbar( 'YDir', 'Reverse' )
title( 'Exact line source solution' )



