% Builds an SSP file for the rotated Munk profile

sspfil = 'MunkRot.ssp';

rmin =    0;
rmax = 5000;
dr   =  200;

zmax = 100000;
dz   =  50000;

r  = [ 0 200 250 ( 400 : dr : rmax ) ] - 1000;
Nr = length( r );

z  = ( 0 : dz : zmax );
Nz = length( z );

% write to an ascii file for BELLHOP
fid = fopen( sspfil, 'w' );
fprintf( fid, '%i \r\n', Nr );
fprintf( fid, '%f  ', r / 1000 ); 
fprintf( fid, '\r\n' );

for iz = 1 : Nz
        for ir = 1 : Nr
            fprintf( fid, '%.2f  ', Munk( r( ir ) + 1000 ) );
        end
        fprintf( fid, '\r\n' );
end
fclose( fid );
   
   
