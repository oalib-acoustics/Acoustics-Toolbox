function [ Nx, Ny, Nz, Segx, Segy, Segz, cmat ] = readssp3d( sspfil )

% Read the 3D SSPFIL used by Bellhop3D
%
% NProf is the number of profiles
% NSSP  is the number of points in depth for a single profile

% open the file

fid = fopen( sspfil, 'r' );

if ( fid == -1 )
   warndlg( 'No SSPFIL file with that name exists', 'Warning' );
end

% read
Nx   = fscanf( fid, '%i', 1 );
Segx = fscanf( fid, '%f', Nx );

Ny   = fscanf( fid, '%i', 1 );
Segy = fscanf( fid, '%f', Ny );

Nz   = fscanf( fid, '%i', 1 );
Segz = fscanf( fid, '%f', Nz );

cmat = zeros( Nz, Ny, Nx );

for iz = 1 : Nz
   cmat2d  = fscanf( fid, '%f', [ Nx, Ny ] );
   cmat( iz, :, : ) = cmat2d';
end

fclose( fid );