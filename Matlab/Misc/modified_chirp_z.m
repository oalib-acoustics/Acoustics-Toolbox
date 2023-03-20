function X = modified_chirp_z( x, M, delta_phi, phi_0 )

% note: the Matlab Signal Processing Toolbox now includes a function czt
% that seems to do the same thing
% The version below, I'm told, has a problem if the vector, x, is too big
% mbp Jan. 2019

% function X = modified_chirp_z( x, M, delta_phi, phi_0 )
% this is a modified chirp z-transform algorithm
% it has been modified to compute the z-plane contour 
% on the unit circle only.
% x is the N point input time series
% the contour begins at |z|=1 and phi = phi_0
% the angular spacing on the unit circle = delta_phi,
% X is returned as an M point complex vector
%
% Reference:
%
% L.R.Rabiner, R.W.Schafer, and C.M.Rader, 
% "The Chirp z-Transform", IEEE Trans. Audio Electroacoustics,
% vol. AU-17, pp 86-92, June 1969.

% J.M.Tattersall NUWC 9/6/96

% mbp: Paul Hursky tells me he thinks this version does not work correctly
% if phi_0 is nonzero.
% Apparently -li * phi_0 should have a sequence (1, ..., N)

% warn the user that the combination of phi_0, delta_phi
% and M has lead to a "wrap" around the unit circle
if phi_0+M*delta_phi > 2*pi
   disp(' Warning:  modified_chirp_z.m')
   disp(' the phase has wrapped around the unit circle' )
end


% the factor 1/2 in the (n^2)/2 term is absorbed 
% into the delta_phi term:
delta_phi = delta_phi/2;

N = length( x );

n_v_sq = ((-(N-1)):1:(M-1)).^2;
v = exp( 1i*delta_phi*n_v_sq);

n_y_sq = (0:N-1).^2;
y = exp(-1i*phi_0)*x.*exp( -1i*delta_phi*n_y_sq);

g = fftfilt( y, v );

g = g(N:N+M-1);

n_X_sq = (0:M-1).^2;

X = g.*exp( -1i*delta_phi*n_X_sq );

return
