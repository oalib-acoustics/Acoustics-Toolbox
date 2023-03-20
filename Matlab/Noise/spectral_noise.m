function NL = spectral_noise( GreenFile, rho_SL_dB, sz, rz, freqIn )
%
% Calculates the NL dues to a uniform sheet of monopoles
%
% This version works off a depth-separated Green's function calculated by SCOOTER
%
% rho_SL_dB source level density (dB expressed per m^2 of the noise sheet)
% GreenFile name of the Green's function file
% sd        source depth (m)
% rd        receiver depth (m) (can be a vector)
% freqIn    source frequency. This is used to pick out one frequency in a
%           broadband run. If the SCOOTER run only had one frequency then
%           this parameter has no effect.
%
% returns
% NL in dB with nrd rows and nrr columns
%
% The source level density assumes source strength reference to the usual
% standard (a source that produces 1 microPascal as 1 m) scaled by the
% density per meter

% The derivation of this noise formula in JKPS assumes the integration is
% done along the real axis. Therefore stabilizing attenuation should not be
% used

% mike porter, 2012

Pos.s.z = sz;   % depth of noise sources
Pos.r.z = rz;   % vector of receiver depths
%Nrd     = length( Pos.r.z );

rho_SL = 10^( rho_SL_dB / 10 ); % source level per m^2 (related to q^2)

%%
% read in the Green's function and interpolate to sd, rd

[ ~, ~, freqVec, ~, ~, PosG, Gtemp ] = read_shd( GreenFile, freqIn );

% Identify the frequency used as the available one closest to the user's request
% If there is only one frequency in the file, then that's the one you'll get
freqdiff = abs( freqVec - freqIn );
[ ~, ifreq ] = min( freqdiff );
freq = freqVec( ifreq );

% calculated k0 used as a scaling factor in Eq. 9.17 JKPS
% Becomes irrelevant because it goes into q2, but q2 is divided by k0 in
% the final noise level
omega = 2 * pi * freq;
c     = 1500;        % nominal sound speed
k0    = omega / c;   % nominal vertical wavenumber

% calculate vector of wavenumbers
k     = 2 * pi * freq ./ PosG.r.r;   % PosG.r.r vector contains phase speeds
atten = 0;   % stabilizing attenuation not allowed !!!

% Get Green's function and interpolate for given source and receiver depth
% interp1 won't go if there is only one point in depth
% FIXME: will fail if there is only one rd as interp1 needs two

if length( PosG.s.z ) > 1
   Gtemp2 = interp1( PosG.s.z, squeeze( Gtemp( 1, :, :, : ) ), Pos.s.z );
   Gtemp2 = squeeze( Gtemp2 );
else
   zs   = Pos.s.z;
   % what is this about? If just one PosG.s.z, need to take it or flag it
   isd  = find( PosG.s.z >= zs );    % index of source depth
   isd  = isd( 1 );
   Gtemp2 = squeeze( Gtemp( 1, isd, :, : ) );
end

% interpolate to get Green's fn G of size Nrd x Nk

G = interp1( PosG.r.z, squeeze( Gtemp2( :, : ) ), Pos.r.z );

%%
% Construct Power-- See Eq. (9.17), p. 667 in JKPS

dk = ( k( end ) - k( 1 ) ) / ( length( k ) - 1 );

%for ird = 1 : Nrd

   % % analytic Green's function
   % % You get the best answer if you set atten = 0
   % gamma = sqrt( k0^2 - ( k - 1i * atten ).^2 );
   % G( ird, : ) = ( exp( 1i * gamma * abs( rd( ird ) - sd ) )  - ...
   %                 exp( 1i * gamma * abs( rd( ird ) + sd ) ) )./ gamma;
   % G2( ird ) = abs( G( ird, : ) ).^2 * ( k - 1i * atten );   % note atten should be 0 for a noise run
%end

G2 = abs( G ).^2 * ( k.' - 1i * atten );   % note atten should be 0 for a noise run
NL = G2 * dk;   % need dk for trapezoidal rule for integral over k

q2 = k0^2 * 4 * pi * rho_SL;          % scale factor to convert a SL density to q2
NL = ( 8 * pi^2 * q2 / k0^2 ) * NL;   % scale factor from Eq. 9.17 in JKPS
NL = NL / ( 4 * pi )^2;               % scale factor to take the SCOOTER definition of g to the jkps one.
% If you simplify the above you get NL = 2 * pi * G2, I think

NL = 10 * log10( NL );   % convert to dB
