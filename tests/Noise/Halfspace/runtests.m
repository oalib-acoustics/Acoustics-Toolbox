% tests for halfspace
% If the source is put at lambda/4 then
% the far field should have a level 13.2 dB higher than rho_SL_dB
% which is 10 log10( 6.6 * pi )
%
% See the Kewley paper for a derivation of that based on a sheet of
% monopoles in a halfspace

% The case is run with stabilizing attenuation turned off so you get the
% 13.2 dB level
% With stabilizing attenuation there's a drop in level with depth

scooter halfspace

%%

GreenFile = 'halfspace.grn';
sd        = 1.25;  % depth of noise sources
rd        = 0 : 1 : 100;
Component = 'P';
rho_SL_dB = 54;   % monopole strength from Kewley for 800 Hz and 40 knots
rho_SL_dB = 0;
freq      = 300;

NL = spectral_noise( GreenFile, rho_SL_dB, sd, rd, freq );

figure
plot( rd, NL, 'LineWidth', 3 )
xlabel( 'Receiver depth (m)' )
ylabel( 'NL (dB)' )
title( 'Halfspace, spectral formula' )

