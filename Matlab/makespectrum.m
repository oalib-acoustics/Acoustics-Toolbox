function [ freq, shat ] = trian( Fmin, Fmax, Nfreq )

% make a triangular spectrum

freq = linspace( Fmin, Fmax, Nfreq );

mid = Nfreq / 2;
midint = floor( mid );

if mid == midint
   shat = [ 0:mid-1, mid-1:-1:0 ];   % Nfreq frequencies
else
   shat = [ 0:mid,   mid-1:-1:0 ];   % Nfreq frequencies
end

for i = 1: size( freq, 2 )
   fprintf(  1, '%6.2f %6.2f 0.0 \n', freq( i ), shat( i ) )
end
