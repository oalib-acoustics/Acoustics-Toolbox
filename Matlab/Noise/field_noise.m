function C = field_noise( ShadeFile, sd, rd, freq )

% Computes a noise covariance matrix by summing columns of the pressure field
% ShadeFile is the name of the shade file
% sd is the vector of source depths
% rd is the vector of receiver depths
%
% Returns:
% C is the nrd x nrd covariance matrix
%
% If you want the noise level in dB use:
%    NL = 10 * log10( diag( C ) ) + rho_SL_dB;

% mike porter Nov. 2017

% read in the pressure array
[ ~, ~, ~, ~, ~, Pos, pressure ] = read_shd( ShadeFile, freq );


%%
% interpolate pressure field to user provided sd, rd

if length( Pos.s.z ) > 1
   Ptemp = interp1( Pos.s.z, squeeze( pressure( 1, :, :, : ) ), sd );
   Ptemp = squeeze( Ptemp );
else
   zs   = sd;

   % find a index of modes.z where the depth is near zs
   isd = find( abs( Pos.s.z - zs ) < eps( max( abs( Pos.s.z ) ) ) );    % index of source depth
   if ( isempty( isd ) )
      errordlg( 'Modes not tabulated at source depth' )
   end
   
   isd  = isd( 1 );
   Ptemp = squeeze( pressure( 1, isd, :, : ) );
end

P = interp1( Pos.r.z, squeeze( Ptemp( :, : ) ), rd );

% loop over range, adding the contribution of each vector into the
% covariance matrix

nrd = length( rd );
C = zeros( nrd, nrd );

for ir = 1 : length( Pos.r.r ) - 1
   
    dr = Pos.r.r( ir + 1 ) - Pos.r.r( ir );
    d = P( :, ir ) * sqrt( 2 * pi * Pos.r.r( ir ) * dr );

    C = C + d * d';
end

