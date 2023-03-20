% Builds a cube with Munk SSP across x-axis constant in y in the z=0 plane
% and polar symmetry elsewhere. 
% In the z=0 plane c(x,y,0) = cMunk(x+1000)for x in [-1000 4000], limited 
% elsewhere.   and c(x,y,z) = c(q,0,0) where q = R*cos(theta), 
% R = sqrt(x^2+y^2+z^2), and theta = atan2(y,x).

sspfil = 'Munk3D.ssp';

Xmin = -10;
Xmax = 5010;
dx   = ( Xmax - Xmin ) / 26;

Ymax = 200000;
dy   = 100000;

Zmax = 3000;
dz   = 1000;

x  = ( Xmin : dx : Xmax );
Nx = length( x );

y  = ( 0 : dy : Ymax );
Ny = length( y );

z  = ( -Zmax : dz : Zmax );
Nz = length( z );


xref = ( Xmin : 10 : Xmax );
cx = Munk( xref );

c = zeros( Nz, Ny, Nx );
X = repmat( x, Ny, 1 );
Y = repmat( y(:), 1, Nx );
theta = atan2( Y, X );
X2Y2 = X.^2 + Y.^2;

for iz = 1 : Nz
    zi = z( iz );
    R = sqrt( zi^2 + X2Y2 );
    q = R .* cos( theta );
    c( iz, :, : ) = interp1( xref, cx, q );
end    
         
c( isnan( c ) ) = 1500;
z = z + Zmax;     % Now z is 0->2*Zmax

% write to an ascii file for BELLHOP3D
fid = fopen( sspfil, 'w' );
fprintf( fid, '%i \r\n', Nx );
fprintf( fid, '%f  ', x / 1000 ); 
fprintf( fid, '\r\n' );

fprintf( fid, '%i \r\n', Ny );
fprintf( fid, '%f  ', y / 1000 );
fprintf( fid, '\r\n' );

fprintf( fid, '%i \r\n', Nz );
fprintf( fid, '%f  ', z );
fprintf( fid, '\r\n' );

for kk = 1 : Nz
    for jj = 1 : Ny
        for ii = 1 : Nx
            fprintf( fid, '%f  ', Munk( x( ii ) ) ); % c( kk, jj, ii ));  % C(z,y,x)
            % fprintf( fid, '%f  ', Munk( z( kk )        ) ); % c( kk, jj, ii ));  % C(z,y,x)

        end
        fprintf( fid, '\r\n' );
    end
end

fclose( fid );
