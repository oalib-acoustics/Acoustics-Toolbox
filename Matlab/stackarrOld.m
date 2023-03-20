
% plot replica time series
% mbp 8/96, Universidade do Algarve

clear all

! bellhop munkts % ssp14_1544

Nrd = 51;
rd = 1:Nrd;	% select the rd to use
deltat = .005;
T = 5.0;	% time window over which impulse response is calculated
Tstart = 19.0;	% start time
c = 1535.0; %1507.5;   % optional reduction velocity

% set up source spectrum
pulsetime = deltat:deltat:T;
s = Ricker( pulsetime, 10 );
shat = fft( s ).';

LoadArrfil
%plot( delay, amp, 'w*' )

% compute the sum

Fmax = 1 / deltat;
nt = T / deltat;	% number of time points

%duration = 0.05; % 0.0025;
%pulsestart = -8 * duration;
%pulseend   = T + pulsestart - deltat;

%timepulse = pulsestart:deltat:pulseend;
%Ntpulse = size( timepulse, 2 );

%svec = s( timepulse, 0.0, duration );
%svech = hilbert( svec );

H = zeros( nt, nr, Nrd );

for ird = rd

for ir = 1:nr
  last = narrmat( ir, rd );

  tstart = Tstart; % rr( ir ) / c - 0.1   % min( delay( ir, 1:last, ird ) )
  tend   = tstart + T - deltat;
  time   = tstart:deltat:tend;

  % compute channel transfer function

  for iarr = 1 : narrmat( ir, ird )

    % compute time index for delay
    it = floor( ( delay( ir, iarr, ird ) - tstart ) / deltat + 1 );
    if ( it >= 1 && it <= nt )
       H( it, ir, ird ) = H( it, ir, ird ) + amp( ir, iarr, ird );
    end

  end   % next arrival, iarr.

end   % next range, ir.

end   % next rd

% form r( t ) as convolution of impulse response with s( t )

rmod = zeros( nt, nr, Nrd );

for ird = rd
  for ir = 1:nr
    ir
    gr = real( H( :, ir, ird ) );
    gi = imag( H( :, ir, ird ) );
    grhat = fft( gr );
    giquad = imag( hilbert( gi ) );
    ghat = grhat + fft( giquad ); 

    % note: ifft below has small im. parts due to round-off
    % - makes sign consistent with KRAKEN/SCOOTER
    rmod( :, ir, ird ) = real( ifft( -ghat .* shat ) );
  end
end

figure; plot( rmod( :, 42 ) )

figure
peak = max( max( rmod( :, 1, : ) ) );
imagesc( time, rd , squeeze( rmod( :, 1, : ) )' )
caxis( [ -peak/5, peak/5 ] ); colorbar
title( 'BELLHOP impulse response' )    

%rds = [ 1 3 4 ];
%for ird = 1:3
%  subplot( 3, 1, ird )
%  plot( abs( squeeze( rmod( 1, :, rds( ird ) ) ) ) )
%end

nrr = nr;
fname = 'model'
eval( [ 'save ' fname ' nrr rr rmod' ] );


  % take log and do some normalization

  rtemp = 20 * log10( r / max( r ) ) + 30;

  % clip
  I = find( rtemp < 0 );
  rtemp( I ) = zeros( size( I ) );
  rmod( ir, :, rd ) = rtemp;
