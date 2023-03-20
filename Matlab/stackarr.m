
% plot replica time series
% mbp 8/96, Universidade do Algarve

clear all

! bellhop ssp16_1500

Fmax = 4000.0;
deltat = 1 / Fmax;

T = 1.0;		% time window over which impulse response is calculated
c = 1507.5 %1535.0;	% optional reduction velocity

% set up source spectrum
pulsetime = deltat : deltat : T;
specfile = 'spta.fft';

% load source spectrum
fid = fopen( specfile );
temp = fscanf( fid, '%f', [ 3, inf ] );

freq = temp( 1, : );
shat = temp( 2, : ) + 1i * temp( 3, : );

shat = [ zeros( 1, 300 ) shat zeros( 1, 3199 ) ]';   % in [ 0, 1000 ] Hz
nfreq = size( shat, 2 );

LoadArrfil

%%
% compute the sum

nt = T / deltat;	% number of time points

for isd = 1:nsd
  H = zeros( nt, nr, nrd );	% clear the transfer function
  for ird = 1:nrd
    for ir = 1:nr
      [ isd, ird, ir ]
      last = narrmat( ir, ird, isd );

      tstart = rr( ir ) / c - 0.1;   % min( delay( ir, 1:last, ird ) )
      tend   = tstart + T - deltat;
      time   = tstart : deltat : tend;

      % compute channel transfer function

      for iarr = 1 : narrmat( ir, ird, isd )

        % compute time index for delay
        Tarr = delay( ir, iarr, ird, isd ) - tstart;
        it = round( Tarr / deltat + 1 );
        % s = Tarr / deltat + 1 - it;		% proportional interval
        if ( it >= 1 && it <= nt )
           H( it,   ir, ird ) = H( it,   ir, ird ) + amp( ir, iarr, ird );
        end

      end   % next arrival, iarr.

    end   % next range, ir.

  end   % next rd
%%
  % form r( t ) as convolution of impulse response with s( t )

  rmod = zeros( nt, nr, nrd );

  for ird = 1 : nrd
    for ir = 1 : nr
      gr     = real( H( :, ir, ird ) );
      gi     = imag( H( :, ir, ird ) );
      grhat  = fft( gr );
      giquad = imag( hilbert( gi ) );
      ghat   = grhat + fft( giquad ); 

      % note: ifft below has small im. parts due to round-off
      % - makes sign consistent with KRAKEN/SCOOTER
      rmod( :, ir, ird ) = real( ifft( -ghat .* shat ) );
    end
  end
%%
  % take log and do some normalization

  rmod = abs( hilbert( rmod ) );

  % normalize

  for ird = nrd
    for ir = 1 : nr
      temp = 20 * log10( rmod( :, ir, ird ) / max( rmod( :, ir, ird ) ) ) + 30;
      I    = find( temp < 0 );
      temp( I ) = zeros( size( I ) );
      rmod( :, ir, ird ) = temp / norm( temp );
    end
  end

  foo   = rmod( 1:4:4000, :, : );
  rmod  = foo;
  nrr   = nr;
  fname = [ 'model' int2str( isd ) ]
  eval( [ 'save ' fname ' nrr rr rmod' ] );
end	% next source depth

%figure
%peak = max( max( rmod2 ) );
%imagesc( time, rd , rmod2 )
%colorbar
%title( 'BELLHOP impulse response' )    

%rds = [ 1 3 4 ];
%for ird = 1:3
%  subplot( 3, 1, ird )
%  plot( abs( squeeze( rmod( 1, :, rds( ird ) ) ) ) )
%end
