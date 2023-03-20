% Run wedge problem to test coupled mode option in KRAKEN
% write out the envfil

freq = 100:1:200;
freq = 25;
nfreq = length( freq );

Nprof = 51;
r = linspace( 0.0, 4000, Nprof );
D = linspace( 200,    0, Nprof );
D( end ) = 1.0;

for ifreq = 1:nfreq
    fid = fopen( 'ENVFIL', 'w' );
    
    for ir = 1:length( r )
        fprintf( fid, [ '''Wedge problem #' int2str( ir ) '''\r\n' ] );
        fprintf( fid, '%8.2f \r\n', freq( ifreq ) );
        fprintf( fid, '2 \r\n' );
        fprintf( fid, '''CVW .'' \r\n' );
        fprintf( fid, '500 0.0 %8.2f \r\n', D( ir ) );
        fprintf( fid, ' 0.0 1500.0 / \r\n' );
        fprintf( fid, ' %8.2f 1500.0 / \r\n', D( ir ) );
        fprintf( fid, '1000 0.0 2000.0 \r\n' );
        fprintf( fid, ' %8.2f 1700.0 0.0 1.5 0.5 0.0 \r\n', D( ir ) );
        fprintf( fid, ' 2000.0 1700.0 / \r\n' );
        fprintf( fid, '''V'', 0.0 \r\n' );
        fprintf( fid, '1400.0  15000.0 \r\n' );
        fprintf( fid, '0.0,			! RMAX (km) \r\n' );
        fprintf( fid, '1		    ! NSD \r\n' );
        fprintf( fid, '0.			! SD(1)  ... \r\n' );
        fprintf( fid, '2001		    ! NRD \r\n' );
        fprintf( fid, '0. 2000.0 /	! RD(1)  ... \r\n' );
    end	% next range
    
    fclose( fid );
    
    % run KRAKEN/FIELD
    
    !krakenc.exe < ENVFIL > foo.prt
    !field.exe < field.flp
    filename = [ 'Wedge' int2str( ifreq ) '.shd' ];
    copyfile( 'SHDFIL', 'wedge.shd' );
    delete MODFIL*
    
end   % next frequency

runplots
delete ENVFIL
