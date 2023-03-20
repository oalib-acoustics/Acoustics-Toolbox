
% Do a broadband model run
% Works by writing a sequence of environment files
% each with a different frequency inserted.
% mbp 1/97

clear all
cd work

nfreq = 101;
Fmin = 30.0;
Fmax = 330.0;

freq = linspace( Fmin, Fmax, nfreq );

for ifreq = 1:nfreq

  % open the environmental file

  filename = 'HST_bb'
  fid = fopen( [ filename '.env' ], 'wt' );

  fprintf( fid, '%s\n', '''dummy title''' );   % title
  fprintf( fid, '%f\n', freq( ifreq ) );   % frequency

  fclose( fid );

  % append the tail to the head and run KRAKEN

  % DOS
  eval( [ '!copy ' filename '.env + ..\HST_K_tail.env > foo.prt|' ] );
  eval( [ '!krakenc.exe <' filename '.env > foo.prt|' ] );
  eval( '!field.exe < ..\field.flp > foo.prt|' );
  fileout = [ 'HST' int2str( ifreq ) ]
  eval( [ '!rename shdfil ' fileout '.shd' ] )


end   % next frequency

cd ..
