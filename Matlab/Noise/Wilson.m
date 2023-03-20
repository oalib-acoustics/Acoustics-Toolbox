function SL = Wilson( f, v, SLtype )
% read wilson.xls and plot the corresponding spectra
% data from Wilson, 1981, Table 1.
% spreadsheet from Chris de Moustier

ss = xlsread( 'wilson.xls' );

% parse the data from the table
freq = ss(      2, 2 : end );    % frequency in Hz
wind = ss( 3 : 15, 1       )';   % wind speed in kn
SLdipole  = ss( 3 : 15, 2 : end );   % spectral density level in dB re 1µPa^2/Hz

SLmonopole = SLdipole - 20 * log10( pi );

switch SLtype
   case 'monopole'
      SL = SLmonopole;
   case 'dipole'
      SL = SLdipole;
end

% figure
% semilogx( freq, SLd, '+-', 'linewidth', 2 )
% grid
% xlabel( 'Frequency (Hz)' )
% ylabel( 'Spectral Density Level (dB re 1 µPa^2/Hz)' )
% title( 'Wilson' )
% axis( [ 10, 1000, 40, 70 ] )

% interpolate to user specified frequency and wind speed

[ x, y ] = meshgrid( log10( freq ), wind );

SL = interp2( x, y, SL, log10( f ), v );
