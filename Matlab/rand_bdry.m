% Generate a random boundary (top or bottom)
% mbp 11/7/00

Dmean = 0.0; % 186.0;	% mean depth

dr		= 1.0;			% sampling in range
Rmax	= 5000.0;		% maximum range (m)
r = dr: dr : Rmax;	% range vector
N = length( r );

dkappa = 1 / Rmax;
kappa_max = 1 / dr;
phase = rand( [ 1  N ] );	% random phases

kappa = dkappa:dkappa:kappa_max;
phi = 5.5e-5 * kappa .^-2.25;	% area A power spectrum from Medwin and Clay, p. 659

A = sqrt( phi ) .* exp( 1i * 2 * pi * phase );		% make amplitude spectrum with random phase

D = Dmean + real( ifft( A ) );						% back to the space domain

figure; plot( r, D );
xlabel( 'Range (m)' )
ylabel( 'Depth (m)' )

% write to file in BELLHOP format (ranges in km, depths in m)

B = [ r'/1000 D' ];

fid = fopen( 'rand.bty', 'w' );
fprintf( fid, '%i \r\n', N );			% number of points
fprintf( fid, '%f %f \r\n', B' );	% range-depth pairs

fclose( fid );


