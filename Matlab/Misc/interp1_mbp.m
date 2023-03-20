function Vq = interp1_mbp( X, V, Xq )

% This is a very slightly modified version of Matlab's interp1 function
% that allows the vector X to have just a single sample point
% The Matlab version requires at least two sample points
%
% When there is just one sample point, then the 'interpolation' returns the
% only value it can, i.e. whatever the single value that is in V

if ( length( X ) == 1 )
    Vq = V;
else
    Vq = interp1( X, V, Xq );
end
