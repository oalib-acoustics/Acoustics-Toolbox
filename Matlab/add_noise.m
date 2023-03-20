function add_noise

% Incorporate the source level and the noise level to an existing set of
% time series
% The receiver timeseries is read from a file and assumed to be based on a
% 0 dB source.
%
% mbp 4/09

clear all

pcmFileRoot = '2211_Sd_100_Arr';   % output file
%pcmFileRoot ='autec_V4kts';   % output file

% read a clean packet for use as a marker
STSFIL = '../sample_band-A_pkt-2_2frames.wav'
[ stsTemp, sample_rate ] = wavread( STSFIL );
sts = stsTemp( 1 : 2 : 5*96000 );  % sub-sample down to 48000/s

% WHOI micromodem, band A
fc    = 10000;
BW    = 10000;
flow  = fc - BW / 2;
fhigh = fc + BW / 2;
T     = 5;   % packets of 5 seconds each
fs    = 48000;
deltat = 1 / fs;

SLdB  = [ 185 185 ];   % source level in dB (total power)
SL    = 10.^( SLdB / 20 );

AdB = linspace( 30, 50, 2 );   % noise amplitude in dB (PSD, not total power)
A   = sqrt( 2 ) * sqrt( fhigh - flow ) * 10.^( AdB / 20 );

% generate the noise time series (filtered gaussian random noise)
% noise_ts is a matrix with each column having the noise vector as a
% different power level
noise_ts = makenoise( fc, BW, T, fs ) * A;

% summary
disp( '*** Noise ***' )
disp( 'Power spectral density = ' )
disp( AdB )

disp( 'Total power =' )
disp( 20 * log10( A ) - 3 )

% check using norm
disp( 'Total power calculated directly from timeseries =' )
disp( 20 * log10( norm( noise_ts ) * sqrt( deltat ) ) - 10 * log10( T ) )

for ird = 1 : 20 % length( Pos.r.z );
   disp( ird )
   eval(  [ 'load ' pcmFileRoot '_Rd_' num2str( ird ) ] );
   for ir = 1 : 100 % length( Pos.r.r );
      rts = rtsmat( :, ir );      % extract the timeseries for this particular receiver
      rts = rts * SL + noise_ts;  % incorporate source level and noise level
      
      for iNL = 1 : length( A )   % loop over noise levels
         
         % write to file
         pcmfile = [ pcmFileRoot '_NL_' num2str( AdB( iNL ) ) '.pcm' ];
         if ( ird == 1 && ir == 1 )
            fid_out = fopen( pcmfile, 'w', 'ieee-le' );
            
            % write marker packet (used to know when modem starts
            % processing)
            fwrite( fid_out, 2^15 * 0.95 * sts / max( abs( sts ) ), 'int16' );
         else
            fid_out = fopen( pcmfile, 'a', 'ieee-le' );
         end
         if ( fid_out == -1 )
            error('Can''t open PCM file for output.');
         end
         
         fwrite( fid_out, 2^15 * 0.95 * rts( :, iNL ) / max( abs( rts( :, iNL ) ) ), 'int16' );
         fclose( fid_out );
      end
      
   end
   
end
