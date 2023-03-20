% plot_noise

Rmax_km = linspace( 1, 100, 50 );   % max range of integration of noise sources

%filename = 'Envs_Like_MBPJOE+HB/000.5_-00.5.mod';
%filename = 'H2OC.mod';
%filename = 'Munk.mod';
%filename = '190.5_040.5.mod'
%filename = 'equator.mod'
ModeFile = 'pekeris.mod';

sd = 1;  % depth of noise sources
%rd = 5;    % vector of receiver depths
%rd =  0: 500: 4500;   % vector of receiver depths (need two values for dphi/dz
rd = 0 : 1 : 100;

SL       = 168.0;   % source level in dB
rho_ship = 3;       % ships per square degree

% following assumes a degree at the equator ...
km_per_degree = 107;   % (is this number right?)
rho_SL = 10^( SL / 10 ) * rho_ship / ( km_per_degree * 1000.0 )^2;  % q2 or source level per m^2

rho_SL_dB = 10 * log10( rho_SL );

% rho_SL_dB = 54;   % monopole strength from Kewley for 800 Hz and 40 knots

Component = 'P';

NL = modal_noise( ModeFile, rho_SL_dB, sd, rd, Rmax_km, Component );

figure
plot( Rmax_km, NL )
xlabel( 'Radius of noise disc (km)' )
ylabel( 'Noise Level (dB)' )


figure
imagesc( Rmax_km, rd, NL )
colorbar
caxis( [ 55 65 ] )
%%
% following is stuff that was used for H2O tests
% Assumes NL is on a linear scale, not dB

% NL( 4, : ) = NL( 2, : ) + 2 * NL( 3, : );   % 'kinetic energy'
% 
% figure
% %plot( Rmax_km, 10 * log10( NL ) )
% range_units = Rmax_km / 4.9;
% 
% 
% NL( 1, : ) = NL( 1, : ) / NL( 1, end );
% NL( 2, : ) = NL( 2, : ) / NL( 2, end );
% NL( 3, : ) = NL( 3, : ) / NL( 3, end );
% NL( 4, : ) = NL( 4, : ) / NL( 4, end );
% 
% plot( range_units, NL )
% xlabel( 'Non-dimensional Radius' )
% 
% title( 'Energy' )
% legend( 'Pressure', 'Vertical * \rho c', 'Horizontal * \rho c', 'Kinetic' )
% % axis( [ 0 200 20 60 ] )
