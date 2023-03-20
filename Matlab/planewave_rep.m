function e = planewave_rep( phone_coords, angles, freq )

% set up matrix, e( angle, phone ) of planewave steering vectors
% usage e = planewave( phone_coords, angles, freq )
% angles are relative to broadside
%
% check whether you want windowing on line 35 ...
%
% mbp, October 99

c     = 1480   % reference sound speed for steering vectors

% make sure angles is a column vector
if ( size( angles, 2 ) ~= 1 )
    angles = angles';
end

% make sure phone_coords is a row vector
if ( size( phone_coords, 1 ) ~= 1 )
    phone_coords = phone_coords';
end

theta_con = 90 - angles;			% convert to angle relative to forward-endfire
theta_rad = angles * pi / 180;	% convert to radians

Nelts = length( phone_coords );

% generate a matrix of steering vectors

omega = 2 * pi * freq;
k0    = omega / c;
e     = exp( 1i * k0 * sin( theta_rad ) * phone_coords );

% window and normalize
% note norm( 1 ) is # of elements for a planewave steering vector

window = hanning( Nelts )';
for itheta = 1:length( angles )
  e( itheta, : ) = e( itheta, : ) / norm( e( itheta, : ) );	% normalize steering vectors
  %e( itheta, : ) = e( itheta, : ) .* window;
end

