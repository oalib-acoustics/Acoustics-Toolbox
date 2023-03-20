% Run tests to verify the point vs. line source option

% point source cases
bellhopM 'freePointB'
plotshd( 'freePointB.shd.mat', 4, 2, 1 );
caxisrev( [ 40 80 ] )

bellhopM 'freePoint_gbtB'
plotshd( 'freePoint_gbtB.shd.mat', 4, 2, 3 );
caxisrev( [ 40 80 ] )

scooterM 'freeSPoint'
plotshd( 'freeSPoint.shd.mat', 4, 2, 5 );
caxisrev( [ 40 80 ] )


%%

% line source cases
bellhopM 'freeLineB'
plotshd( 'freeLineB.shd.mat', 4, 2, 2 );
caxisrev( [ 0 25 ] )

bellhopM 'freeLine_gbtB'
plotshd( 'freeLine_gbtB.shd.mat', 4, 2, 4 );
caxisrev( [ 0 25 ] )

scooterM 'freeSLine'
plotshd( 'freeSLine.shd.mat', 4, 2, 6 );
caxisrev( [ 0 25 ] )

%%

% exact solution

freq = 5;
k0 = 2 * pi * freq / 1500.0;
zs = 3000;
zr = 0:25:5000;
rr = 25:25:10000;

ppoint = zeros( length( zr ), length( rr ) );
pline  = zeros( length( zr ), length( rr ) );

for iz = 1: length( zr )
   rmat1 = sqrt( rr.^2 + (zr( iz ) - zs )^2 );
   rmat2 = sqrt( rr.^2 + (zr( iz ) + zs )^2 );
   ppoint( iz, : ) = exp( 1i * k0 * rmat1 ) ./ rmat1;
   pline(  iz, : ) = besselh( 0, k0 * rmat1 );
end

pline = sqrt( pi / 2 ) * pline; %normalization to roughly match acoustics toolbox

subplot( 4, 2, 7 )
pcolor( rr, zr, -20 * log10( abs( ppoint ) ) ); shading interp; view( 0, -90 )
caxisrev( [ 40 80 ] )
title( 'Exact point source solution' )

subplot( 4, 2, 8 )
pcolor( rr, zr, -20 * log10( abs( pline ) ) ); shading interp; view( 0, -90 )
caxisrev( [  0 25 ] )
title( 'Exact line source solution' )

