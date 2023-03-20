
% compares noise using various approaches
% BELLHOP and SCOOTER are using a logarithmic spacing in range
% with about 1000 receivers
%
% KRAKEN and KRAKENC are using a uniform spacing in range with 10,000
% receivers
% The receiver range spacing is not used for the analytic formulas
% BELLHOP as usual has problems for receivers very close to the boundaries
%
% KRAKENC is the only model that produces a TL(z) that is not slightly
% noisy when extracting it from the pressure field (as opposed to using the
% analytic formula)
% It's not clear why this is. SCOOTER uses the same sampling as KRAKENC.
% Tests with the phase speed limits showed no sensitivity.
% A run with the full Hankel transform also produce the same noise in the
% result.
%
% The BELLHOP noise level is too high for deeper depths. It looks like
% there is a string of caustics at depth, which may explain the high levels

global units
units = 'km';

freq      = 300;
v         = 10;   % wind speed in knots
rho_SL_dB = Kewley( freq, v, 'monopole' );   % monopole strength due to wind

sd        = 1.25;  % depth of noise sources (lambda/4)
freq      = 300;   % frequency
%%
% BELLHOP run

bellhop aetB_TL

% figure; plotshd aetB_TL.shd
% caxis( [ 80 130 ] )
% print -dpng BELLHOP_TL

%%
% noise from shd

bellhop aetB
rd = 0 : 10 : 6000;
C  = field_noise( 'aetB.shd', sd, rd, freq );
NLB = 10 * log10( diag( C ) ) + rho_SL_dB;

figure
plot( NLB, rd, 'LineWidth', 3 );
grid
set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
title( 'BELLHOP Noise Level' )
xlabel( 'dB' )
ylabel( 'Depth (m)' )
axis( [ 50 80 0 6000 ] )
drawnow
print -dpng BELLHOP_NL

%%
% KRAKEN run

kraken aetK_TL

% figure; plotshd aetK_TL.shd.mat
% caxis( [ 80 130 ] )
% print -dpng KRAKEN_TL

% noise from shd

kraken aetK
rd = 0 : 10 : 6000;
C  = field_noise( 'aetK.shd.mat', sd, rd, freq );
NLK = 10 * log10( diag( C ) ) + rho_SL_dB;

figure
plot( NLK, rd, 'LineWidth', 3 );
grid
set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
title( 'KRAKEN Noise Level' )
xlabel( 'dB' )
ylabel( 'Depth (m)' )
axis( [ 50 80 0 6000 ] )
drawnow
print -dpng KRAKEN_NL
%%
% KRAKENC run

krakenc aetC_TL

% figure; plotshd aetC_TL.shd.mat
% caxis( [ 80 130 ] )
% print -dpng KRAKENC_TL

% noise from shd

krakenc aetC
rd = 0 : 10 : 6000;
C  = field_noise( 'aetC.shd.mat', sd, rd, freq );
NLC = 10 * log10( diag( C ) ) + rho_SL_dB;

figure
plot( NLC, rd, 'LineWidth', 3 );
grid
set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
title( 'KRAKENC Noise Level' )
xlabel( 'dB' )
ylabel( 'Depth (m)' )
axis( [ 50 80 0 6000 ] )
drawnow
print -dpng KRAKENC_NL

%%
% SCOOTER run

scooter aetS_TL

% figure; plotshd aetS_TL.shd.mat
% caxis( [ 80 130 ] )
% print -dpng SCOOTER_TL

% noise from shd

scooter aetS
rd = 0 : 10 : 6000;
C  = field_noise( 'aetS.shd.mat', sd, rd, freq );
NLS = 10 * log10( diag( C ) ) + rho_SL_dB;

figure
plot( NLS, rd, 'LineWidth', 3 );
grid
set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
title( 'SCOOTER Noise Level' )
xlabel( 'dB' )
ylabel( 'Depth (m)' )
axis( [ 50 80 0 6000 ] )

print -dpng SCOOTER_NL

%%
% overplot results from all the models as extracted from the SHDFIL

plotshd( 'aetB_TL.shd', 2, 2, 1 )
caxis( [ 80 130 ] )

plotshd( 'aetK_TL.shd.mat', 2, 2, 2 )
caxis( [ 80 130 ] )

plotshd( 'aetC_TL.shd.mat', 2, 2, 3 )
caxis( [ 80 130 ] )

plotshd( 'aetS_TL.shd.mat', 2, 2, 4 )
caxis( [ 80 130 ] )

figure
plot( NLB, rd, 'k', 'LineWidth', 3 ); hold on
plot( NLK, rd, 'b', 'LineWidth', 3 ); hold on
plot( NLC, rd, 'r', 'LineWidth', 3 ); hold on
plot( NLS, rd, 'g', 'LineWidth', 3 ); hold on

grid
set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
title( 'Noise Level' )
xlabel( 'dB' )
ylabel( 'Depth (m)' )
axis( [ 50 80 0 6000 ] )

legend( 'Bellhop', 'Kraken', 'Krakenc', 'Scooter' )


%%
% Spectral noise vs. depth (analytic)
% Note that SCOOTER has to be run with stabilizing attenuation disabled
% ( TopOpt( 7 : 7 ) = '0' )

GreenFile = 'aetSNoiseRun.grn';
scooter_nofield aetSNoiseRun

freqIn = 0;   % frequency (irrelevant if there is only one frequency in the SCOOTER run
rd = 0 : 10 : 6000;
NL = spectral_noise( GreenFile, rho_SL_dB, sd, rd, freqIn );

figure
plot( NL, rd, 'r', 'LineWidth', 3 );
grid
set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
title( 'SCOOTER Noise Level (analytic)' )
xlabel( 'dB' )
ylabel( 'Depth (m)' )
axis( [ 50 80 0 6000 ] )

print -dpng SCOOTER_NL_spectral

%%
% KRAKEN Modal noise vs. depth using full matrix
% modal_noise_full takes a really long time to run ...

%Rmax_km = linspace( 1, 1e4, 2 );   % max range of integration of noise sources
Rmax_km = 1e9;   % max range of integration of noise sources

ModeFile  = 'aetK.mod';
ModeFileC = 'aetC.mod';

Component = 'P';
rd = 0 : 20 : 6000;

NL  = modal_noise_full( ModeFile,  rho_SL_dB, sd, rd, freq, Rmax_km, Component );
NLC = modal_noise_full( ModeFileC, rho_SL_dB, sd, rd, freq, Rmax_km, Component );

% we just take the value from the matrix NL for the largest range

figure
plot( NL(  :, end ), rd, 'g', 'LineWidth', 3 ); hold on
plot( NLC( :, end ), rd, 'b', 'LineWidth', 3 );

ylabel( 'Receiver depth (m)' )
xlabel( 'NL (dB)' )
grid
set(gca,'YDir','reverse')
legend( 'Kraken', 'Krakenc' )
title( 'Kraken Noise Level (analytic)' );
xlabel( 'dB' )
ylabel( 'Depth (m)' )
axis( [ 50 80 0 6000 ] )

print -dpng KRAKEN_NL_modal_full

%%
% KRAKEN Modal noise vs. depth using diagonal terms only

Rmax_km = 1e9;   % max range of integration of noise sources

ModeFile  = 'aetK.mod';
ModeFileC = 'aetC.mod';

Component = 'P';
rd = 0 : 10 : 6000;

NL  = modal_noise_diag( ModeFile,  rho_SL_dB, sd, rd, freq, Rmax_km, Component );
NLC = modal_noise_diag( ModeFileC, rho_SL_dB, sd, rd, freq, Rmax_km, Component );

figure
plot( NL(  :, end ), rd, 'g', 'LineWidth', 3 ); hold on
plot( NLC( :, end ), rd, 'b', 'LineWidth', 3 );

ylabel( 'Receiver depth (m)' )
xlabel( 'NL (dB)' )
grid
set(gca,'YDir','reverse')
legend( 'Kraken', 'Krakenc' )
title( 'Kraken/Krakenc Noise Level (analytic, diag. only)' );
xlabel( 'dB' )
ylabel( 'Depth (m)' )
axis( [ 50 80 0 6000 ] )

print -dpng KRAKEN_NL_modal_diag

%%
% Noise directionality
% beamformed look at the noise field
% The original version used a 32-element array
% This version has higher sampling

%rd = linspace( 802.5, 880, 32 );   % 32-element array, 2.5 m spacing
rd = linspace( 800, 880, 81 );   % 32-element array, 2.5 m spacing

bellhop aet_VLA_B
C_B  = field_noise( 'aet_VLA_B.shd', sd, rd, freq );
C_B = C_B * 10^( rho_SL_dB / 10 );   % scale up by the noise level

kraken aet_VLA_K
C_K  = field_noise( 'aet_VLA_K.shd.mat', sd, rd, freq );
C_K = C_K * 10^( rho_SL_dB / 10 );   % scale up by the noise level

krakenc aet_VLA_C
C_C  = field_noise( 'aet_VLA_C.shd.mat', sd, rd, freq );
C_C = C_C * 10^( rho_SL_dB / 10 );   % scale up by the noise level

scooter aet_VLA_S
C_S  = field_noise( 'aet_VLA_S.shd.mat', sd, rd, freq );
C_S = C_S * 10^( rho_SL_dB / 10 );   % scale up by the noise level

%%

% To get the result presented at the NRL Noise Workshop you need to have
% the Hanning window turned on in planewave_rep !!!

angles = -90 : 1 : 90;
e = planewave_rep( rd, angles, freq );   % construct steering vectors
e = e.';      % e is now a matrix of nrd x nangles

nangles = length( angles );

Beam_pattern = zeros( nangles, nangles );

NoiseB = zeros( 1, nangles );
NoiseK = zeros( 1, nangles );
NoiseC = zeros( 1, nangles );
NoiseS = zeros( 1, nangles );

for iangle = 1 : nangles
   norm_e = norm( e( :, iangle ) );
   R = e( :, iangle ) * e( :, iangle )' / norm_e;   % fixed look direction
   % Beam_pattern is real, symmetric, but we'll compute both upper and lower parts
   for iangle2 = 1 : nangles   
      norm_e2 = norm( e( :, iangle2 ) );
      Beam_pattern( iangle, iangle2 ) = real( e( :, iangle2 )' * R * e( :, iangle2 ) ) / norm_e2;      % N = e' R e
   end
   
  NoiseB( iangle ) = 10 * log10 ( real( e( :, iangle )' * C_B * e( :, iangle ) ) / norm_e^2 );      % N = e' C e
  NoiseK( iangle ) = 10 * log10 ( real( e( :, iangle )' * C_K * e( :, iangle ) ) / norm_e^2 );      % N = e' C e
  NoiseC( iangle ) = 10 * log10 ( real( e( :, iangle )' * C_C * e( :, iangle ) ) / norm_e^2 );      % N = e' C e
  NoiseS( iangle ) = 10 * log10 ( real( e( :, iangle )' * C_S * e( :, iangle ) ) / norm_e^2 );      % N = e' C e
  %NoiseS( iangle ) = 10 * log10 ( 1 ./ real( e( :, iangle )' * inv( C_S ) * e( :, iangle ) ) * norm_e );      % MVDR
end
disp( 'should we normalize by norm( e )?' )

%ii = find( Beam_pattern < 1e-100 );
%Beam_pattern( ii ) = 1e-100;

Beam_pattern_dB = 10 * log10( Beam_pattern );
   
figure
plot( angles, NoiseB, 'k', 'LineWidth', 3 ); hold on
plot( angles, NoiseK, 'b', 'LineWidth', 3 ); hold on
plot( angles, NoiseC, 'r', 'LineWidth', 3 ); hold on
plot( angles, NoiseS, 'g', 'LineWidth', 3 ); hold on

legend( 'Bellhop', 'Kraken', 'Krakenc', 'Scooter' )
xlabel( 'conical angle' )
ylabel( 'NL (dB)' )
title( 'Noise Level vs. angle' )

print -dpng NLvsAngle

%%
% Plot noise directionality by converting beam power

%Beam_pattern = eye( size( Beam_pattern ) );

% NoiseSteradian    = real( Beam_pattern \ ( 10 .^ ( real( NoiseB )/ 10 )' ) );
% NoiseSteradian_dB = real( 10 * log10( NoiseSteradian ) );
% 
% figure
% plot( angles, NoiseSteradian_dB )
% 
% 


