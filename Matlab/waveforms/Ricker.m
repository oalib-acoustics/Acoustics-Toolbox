function s = Ricker( time, F )

% make a Ricker wavelet
% peak at F, support [0, 2F]
% F = nominal source frequency

u = 2 * pi * F * time - 8;   % dimensionless time
s = 0.5 * ( 0.25 * u.^2 - .5 ) * sqrt( pi ) .* exp( -0.25 * u.^2 );

deltaf = 1 / max( time );
deltat = time( 2 ) - time( 1 ); 
Fmax = 1 / deltat;
freq = 0:deltaf:Fmax-deltaf;
shat = fft( s );

for i = 2: 151
%   fprintf(  1, '%6.2f %6.2f %6.2f \n', freq( i ), real( shat( i ) ), imag( shat( i ) ) )
end
