
% compares noise from the modal and spectral integral formulas

% run KRAKEN; done this way since we only want the modes, not the file

figure
sd = 0.5;  % depth of noise sources
  rd = 0:1:50;
rho_SL_dB = 0;

for ifile = 1 : 4 % 1 : 9
  filelist = ...
  [  'HS_100Hz      '
     'HS_200Hz      '
     'HS_400Hz      '
     'HS_800Hz      '
      ];
  filename = deblank( filelist( ifile,: ) );


  %%
  % Spectral noise vs. depth
  % Note that SCOOTER has to be run with stabilizing attenuation disabled
  % ( TopOpt( 6 : 6 ) = '0' )

  ScootFile = [ filename 'S' ];
  GreenFile = [ filename 'S.grn' ];
  
%   runscooter = which( 'scooter.exe' );
% 
%   if ( isempty( runscooter ) )
%     error( 'scooter.exe not found in your Matlab path' )
%   else
%     eval( [ '! "' runscooter '" ' ScootFile ] );
%   end

ScootFile
  scooter( ScootFile )
  NL = spectral_noise( GreenFile, rho_SL_dB, sd, rd );
  
  pp = plot( NL, rd );
  set(pp,'LineWidth',3);
  ll = legend( 'Spectral' );
  % ll = legend( 'NM Full', 'NMC Full', 'NM Diag', 'NMC Diag' );
  set(ll,'Location','southeast')
  grid
  hold on
  drawnow 
  
  %%
  NL = noise_shd( [ filename 'S' '.shd.mat' ], rho_SL_dB, sd, rd );
  
  pp = plot( NL, rd );
  set(pp,'LineWidth',3);
  ll = legend( 'Spectral' );
  % ll = legend( 'NM Full', 'NMC Full', 'NM Diag', 'NMC Diag' );
  set(ll,'Location','southeast')
  grid
  hold on
  drawnow 
end


