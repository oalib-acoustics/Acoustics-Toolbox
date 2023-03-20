function SL = Kewley( f, v, SLtype )

% Kewley, et al. curves for noise source level density
% usage:
%    SL = Kewley( f, v, SLtype )
%    f is the frequency  (Hz)    (can be a vector)
%    v is the wind speed (knots) (can be a vector)
%    SLtype is 'monopole' or 'dipole'
%

% mike porter, March 2013

% data is read off the plot in Fig. 6 ('vertical SL')
% We subtracted 6 dB to convert it to a monopole strength for a sheet at lambda/4

freq = [ 35, 600, 1500 ];   % frequency (Hz)
wind = [ 5, 10, 20, 30, 40 ];   % wind speed in knots

% note SL matrix is transposed
% rows of SL correspond to curves at a single wind speed
SLvert = [ ...
   47, 50, 55, 59, 63;    %   35 Hz
   38, 45, 55, 59, 63;    %  600 Hz
   36, 41, 47, 51, 54 ]'; % 1500 Hz

SLmonopole = SLvert - 6;
SLdipole   = SLmonopole + 20 * log10( pi );

switch SLtype
   case 'monopole'
      SL = SLmonopole;
   case 'dipole'
      SL = SLdipole;
end

[ x, y ] = meshgrid( log10( freq ), wind );

SL = interp2( x, y, SL, log10( f ), v );
