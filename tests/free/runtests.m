% Run tests to verify the point vs. line source options
% For the paraxial beam tracing, I had been using a small value for RLoop
% Then, the field looked perfect in the first few km and dropped off to
% quickly in the far field.
% RLoop has been increased to 5 km and the results look perfect
% Increasing RLoop means the beams are wider in the near field, but diverge
% more slowly, producing the minimum width at the 5 km range.

% You need to run the Paraxial beam cases with a ray-centered beam to get a
% good representation of the field for the steep angles near the origin


%%
% BELLHOP point source cases
bellhop 'freePointB'
plotshd( 'freePointB.shd', 3, 2, 1 );
caxisrev( [ 60 80 ] )

bellhop 'freePoint_gbtB'
plotshd( 'freePoint_gbtB.shd', 3, 2, 3 );
caxisrev( [ 60 80 ] )

bellhop 'freePoint_ParaxialB'
plotshd( 'freePoint_ParaxialB.shd', 3, 2, 5 );
caxisrev( [ 60 80 ] )

%%

% BELLHOP line source cases
bellhop 'freeLineB'
plotshd( 'freeLineB.shd', 3, 2, 2 );
caxisrev( [ 10 25 ] )

bellhop 'freeLine_gbtB'
plotshd( 'freeLine_gbtB.shd', 3, 2, 4 );
caxisrev( [ 10 25 ] )

bellhop 'freeLine_ParaxialB'
plotshd( 'freeLine_ParaxialB.shd', 3, 2, 6 );
caxisrev( [ 10 25 ] )

%%

% SCOOTER

scooter 'freeSPoint'
plotshd( 'freeSPoint.shd.mat', 2, 2, 1 );
caxisrev( [ 60 80 ] )

scooter 'freeSLine'
plotshd( 'freeSLine.shd.mat', 2, 2, 2 );
caxisrev( [ 10 25 ] )

%%
% exact (analytic) solution

freq = 5;
k0 = 2 * pi * freq / 1500.0;
zs = 3000;
zr = 0:25:5000;
rr = 25:25:10000;
for iz = 1: length( zr )
   rmat1 = sqrt( rr.^2 + (zr( iz ) - zs )^2 );
   rmat2 = sqrt( rr.^2 + (zr( iz ) + zs )^2 );
   ppoint( iz, : ) = exp( 1i * k0 * rmat1 ) ./ rmat1;
   pline(  iz, : ) = besselh( 0, k0 * rmat1 );
end

pline = sqrt( pi / 2 ) * pline; %normalization to roughly match acoustics toolbox

subplot( 2, 2, 3 )
pcolor( rr / 1000, zr, -20 * log10( abs( ppoint ) ) ); shading interp; view( 0, -90 )
caxisrev( [ 60 80 ] )
title( 'Exact point source solution' )
xlabel( 'Range (km)' )
ylabel( 'Depth (m)' )

subplot( 2, 2, 4 )
pcolor( rr / 1000, zr, -20 * log10( abs( pline ) ) ); shading interp; view( 0, -90 )
caxisrev( [  10 25 ] )
title( 'Exact line source solution' )
xlabel( 'Range (km)' )
ylabel( 'Depth (m)' )

