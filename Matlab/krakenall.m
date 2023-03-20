
% run kraken for every environmental file in a folder

runkraken = which( 'kraken.exe' );

DirInfo = dir( '*.env' );
Nfiles = length( DirInfo );

h = waitbar( 1 / Nfiles, 'Please wait ...' );

for ifile = 1 : Nfiles
    [ ~, filename, ~ ] = fileparts( DirInfo( ifile ).name );
    disp( filename )
    waitbar( ifile / Nfiles );
    eval( [ '! "' runkraken '" ' filename ] );
end

close( h )   % delete the waitbar