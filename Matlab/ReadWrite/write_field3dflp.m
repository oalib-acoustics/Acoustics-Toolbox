
envfil = 'KoreanSea_3D';
btyfil = 'KoreanSea.bty';
flpfil = [ envfil '.flp' ];
model = 'BELLHOP3D';

[ TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, ~ ] = read_env( envfil, model );

Sx  = Pos.s.x;
Sy  = Pos.s.y;
Nsx = Pos.Nsx;
Nsy = Pos.Nsy;

% We already made all the X_Y.mod files with Make_1Deg_Envs.m.
% d = dir('*Bathy.mat'); mn=char(d.name); xyzfile = mn; %bathymetry file
% load( xyzfile ); %bathymetry file
%
% nx = length( Bathy.Lon );
% ny = length( Bathy.Lat );

% set up the field3d.flp file

fid = fopen( flpfil, 'wt' );

fprintf( fid, '''%s'' ! Title', TitleEnv );
fprintf( fid, '''STDFM''	      	! OPT' );
% fprintf( fid, '''%4s'' ! Option', Option );

fprintf( fid, '999999   ! Mlimit (number of modes to include)' );

%%

% source x-coordinates
fprintf( fid,  '%d                 ! Nsx' , Nsx );
fprintf( fid,  '%d %d          /   ! Sx( 1 : Nsx )  coordinates of source (km)', ...
    Sx( 1 ), Sx( end ) );

% source y-coordinates
fprintf( fid,  '%d                 ! Nsy' , Nsy );
fprintf( fid,  '%d %d          /   ! Sy( 1 : Nsy )  coordinates of source (km)', ...
    Sy( 1 ), Sy( end ) );

% source depths

fprintf( fid, '%5i \t \t \t \t ! NSD', length( Pos.s.z ) );

if ( length( Pos.s.z ) > 2 && equally_spaced( Pos.s.z ) )
    fprintf( fid, '    %6f  ', Pos.s.z( 1 ), Pos.s.z( end ) );
else
    fprintf( fid, '    %6f  ', Pos.s.z );
end

fprintf( fid, '/ \t ! SD(1)  ... (m)' );

% receiver depths

fprintf( fid, '%5i \t \t \t \t ! NRD', length( Pos.r.z ) );

if ( length( Pos.r.z ) > 2 && equally_spaced( Pos.r.z ) )
    fprintf( fid, '    %6f  ', Pos.r.z( 1 ), Pos.r.z( end ) );
else
    fprintf( fid, '    %6f  ', Pos.r.z );
end

fprintf( fid, '/ \t ! RD(1)  ... (m)' );

% receiver ranges

fprintf( fid, '%5i \t \t \t \t ! NRR', length( Pos.r.r ) );

if ( length( Pos.r.r ) > 2 && equally_spaced( Pos.r.r ) )
    fprintf( fid, '    %6f  ', Pos.r.r( 1 ), Pos.r.r( end ) );
else
    fprintf( fid, '    %6f  ', Pos.r.r );
end
fprintf( fid, '/ \t ! RR(1)  ... (km)' );

% receiver bearings
fprintf( fid, '37              ' );
fprintf( fid, '0.0 360.0 /        ! NTHETA THETA(1:NTHETA) (degrees)' );

%%
% Read the bathymetry
% note that Bathy.depth has been transposed
[ Bathy.X, Bathy.Y, Bathy.depth, NbtyPtsx, NbtyPtsy ] = readbty3d( btyfil );

% subsample
iskip = 10;
Bathy.X     = Bathy.X(     1 : iskip : NbtyPtsx );
Bathy.Y     = Bathy.Y(     1 : iskip : NbtyPtsy );
Bathy.depth = Bathy.depth( 1 : iskip : NbtyPtsy, 1 : iskip : NbtyPtsx );

NbtyPtsx = length( Bathy.X );
NbtyPtsy = length( Bathy.Y );

%%
% nodes

nnodes = NbtyPtsx * NbtyPtsy;
fprintf( fid, '%5i', nnodes );

% The order of these do-loops must match the order of the do-loops for the nodes.
for iy = 1 : NbtyPtsy
    for ix = 1 : NbtyPtsx
        if ( Bathy.depth( iy, ix ) > 0 )
            modfil = [ '''StB_' num2str( Bathy.X( ix ),'%07.1f' ) '_' num2str( Bathy.Y( iy ),'%07.1f' ) '''' ];
            
        else
            modfil = ' ''DUMMY'' ';
        end
        fprintf( fid, '%8.2f %8.2f %s', Bathy.X( ix ), Bathy.Y( iy ), modfil );
    end
end

% elements
% done in triangle pairs

nelts = 2 * ( NbtyPtsx - 1 ) * ( NbtyPtsy - 1 );
fprintf( fid, '%5i', nelts );

inode = 1;

% The order of these do-loops must match the order of the do-loops for
% finding the mod files.
for iy = 1 : NbtyPtsy - 1
    for ix = 1 : NbtyPtsx - 1
        fprintf( fid, '%5i %5i %5i', inode,     inode + 1,  inode + NbtyPtsx     );
        fprintf( fid, '%5i %5i %5i', inode + 1, inode + NbtyPtsx, inode + NbtyPtsx + 1 );
        inode = inode + 1;
    end
    inode = inode + 1;
end

fclose( fid );
