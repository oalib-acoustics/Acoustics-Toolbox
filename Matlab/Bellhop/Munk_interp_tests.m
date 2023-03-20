

r = linspace( 0, 100000, 501 )';
z = linspace( 0,  5000,  501 )';
dr = r( 2 ) - r( 1 )
dz = z( 2 ) - z( 1 )

c0   = 1500.0;
coef = 0.00737;

DsDr = -500/100000;
DsDr = 0
'range independent !!!'

% pre-allocate
c = zeros( length( z ), length( r ) );
cz = c;
cr = c;
czz = c;
crr = c;
crz = c;

% analytic derivatives

for ir = 1:length( r )
   sspaxis = 1300.0 + DsDr * r( ir );
   
   xt     =  2 * ( z - sspaxis ) / sspaxis;
   DxDz   =  2 / sspaxis;
   DxDr   = -2 * z / sspaxis^2 * DsDr;
   DxDrDz = -2 / sspaxis^2 * DsDr;
   DxDrDr =  4 * z / sspaxis^3 * DsDr^2;
   
   c(   :, ir )  = c0 * ( 1.0 + coef * ( xt - 1.0 + exp( -xt ) ) );
   cz(  :, ir )  = c0 * coef * ( 1 - exp( -xt ) ) .* DxDz;
   cr(  :, ir )  = c0 * coef * ( 1 - exp( -xt ) ) .* DxDr;
   czz( :, ir )  = c0 * coef * exp( -xt ) .* DxDz .^2;
   crz( :, ir )  = c0 * coef * ( exp( -xt ) .* DxDz .* DxDr + ( 1 - exp( -xt ) ) .* DxDrDz );
   crr( :, ir )  = c0 * coef * ( exp( -xt ) .* DxDr .^2     + ( 1 - exp( -xt ) ) .* DxDrDr );
end

figure;
subplot( 2, 3, 1 ); imagesc( c   ); colorbar
subplot( 2, 3, 2 ); imagesc( cz  ); colorbar
subplot( 2, 3, 3 ); imagesc( cr  ); colorbar
subplot( 2, 3, 4 ); imagesc( crr ); colorbar
subplot( 2, 3, 5 ); imagesc( crz ); colorbar
subplot( 2, 3, 6 ); imagesc( czz ); colorbar

% numerical derivatives

[ cr,  cz  ] = gradient( c,  dr, dz );
[ crr, crz ] = gradient( cr, dr, dz );
[ czr, czz ] = gradient( cz, dr, dz );

figure;
subplot( 2, 3, 1 ); imagesc( c   ); colorbar
subplot( 2, 3, 2 ); imagesc( cz  ); colorbar
subplot( 2, 3, 3 ); imagesc( cr  ); colorbar
subplot( 2, 3, 4 ); imagesc( crr ); colorbar
subplot( 2, 3, 5 ); imagesc( crz ); colorbar
subplot( 2, 3, 6 ); imagesc( czz ); colorbar

save sspgrid r z c cz cr crr crz czz
