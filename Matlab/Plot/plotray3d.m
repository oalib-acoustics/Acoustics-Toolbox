function varargout = plotray3d( rayfil )

% Plot the RAYfil produced by Bellhop3D
% usage: plotray( rayfil )
% where rayfil is the ray file (extension is optional)
% e.g. plotray( 'foofoo' )
%
% MBP July 1999

global units jkpsflag

if ( strcmp( rayfil, 'RAYFIL' ) == 0 && ~contains( rayfil, '.ray' ) )
    rayfil = [ rayfil '.ray' ]; % append extension
end

% plots a BELLHOP ray file

fid = fopen( rayfil, 'r' );   % open the file
if ( fid == -1 )
    disp( rayfil );
    error( 'No ray file exists; you must run BELLHOP first (with ray ouput selected)' );
end

% read header stuff

TITLE       = fgetl(  fid );
FREQ        = fscanf( fid, '%f', 1 );
Nsxyz       = fscanf( fid, '%f', 3 );
NBeamAngles = fscanf( fid, '%i', 2 );

DEPTHT      = fscanf( fid, '%f', 1 );
DEPTHB      = fscanf( fid, '%f', 1 );
Type        = fgetl( fid );
Type        = fgetl( fid );

Nsx    = Nsxyz( 1 );
Nsy    = Nsxyz( 2 );
Nsz    = Nsxyz( 3 );

Nalpha = NBeamAngles( 1 );
Nbeta  = NBeamAngles( 2 );

% Extract letters between the quotes
nchars = strfind( TITLE, '''' );   % find quotes
TITLE  = [ TITLE( nchars( 1 ) + 1 : nchars( 2 ) - 1 ) blanks( 7 - ( nchars( 2 ) - nchars( 1 ) ) ) ];
TITLE  = deblank( TITLE );  % remove whitespace

% Extract letters between the quotes
nchars = strfind( Type, '''' );   % find quotes
Type   = Type( nchars( 1 ) + 1 : nchars( 2 ) - 1 );

% read rays

% set up axis lengths for publication
if ( jkpsflag )
    set( gcf, 'Units', 'centimeters' )
    set( gca, 'ActivePositionProperty', 'Position', 'Units', 'centimeters' )
    
    set( gca, 'Position', [ 2 2 14.0  7.0 ] )
    %set( gcf, 'PaperPosition', [ 3 3 19.0 10.0 ] )
end

% set( gca, 'YDir', 'Reverse' )   % plot with depth-axis positive down

xlabel( 'Range, x (m)' )
ylabel( 'Range, y (m)' )
if ( strcmp( units, 'km' ) )
    xlabel( 'Range, x (km)' )
    ylabel( 'Range, y (km)' )
end

zlabel( 'Depth (m)' )
title( TITLE )
hold on

% axis limits
xmin = +1e9;
xmax = -1e9;
ymin = +1e9;
ymax = -1e9;
zmin = +1e9;
zmax = -1e9;

for isx = 1 : Nsx
    for isy = 1 : Nsy
        %[ isx, isy ]
        for ibeam2 = 1 : Nbeta
            for ibeam = 1 : Nalpha
                
                alpha0    = fscanf( fid, '%f', 1 );
                nsteps    = fscanf( fid, '%i', 1 );
                NumTopBnc = fscanf( fid, '%i', 1 );
                NumBotBnc = fscanf( fid, '%i', 1 );
                if isempty( nsteps ); break; end
                
                switch Type( 1 : 2 )
                    case 'rz'
                        ray = fscanf( fid, '%f', [2 nsteps] );
                        r = ray( 1, : );
                        z = ray( 2, : );
                        
                        y = -r;
                        x = zeros( size( r ) );
                    case 'xy'
                        ray = fscanf( fid, '%f', [3 nsteps] );
                        x =  ray( 1, : );
                        y =  ray( 2, : );
                        z =  ray( 3, : );
                end
                
                if ( strcmp( units, 'km' ) )
                    x = x / 1000;   % convert to km
                    y = y / 1000;
                end
                
                %lincol = 'kbgrcmy';
                %ii = NumBotBnc;
                %ii = mod( ii, 3 ) + 1;
                %plot( r, z, lincol(ii) );
                if NumTopBnc >= 1 && NumBotBnc >= 1
                    plot3( x, y, z, 'k' )    % hits both boundaries
                elseif NumBotBnc >= 1
                    plot3( x, y, z, 'b'  )	% hits bottom only
                elseif NumTopBnc >= 1
                    plot3( x, y, z, 'g'  )	% hits surface only
                else
                    plot3( x, y, z, 'r' )
                end
                
                %                 switch ibeam2
                %                     case 1
                %                         plot3( x, y, z, 'k' )    % hits both boundaries
                %                     case 2
                %                         plot3( x, y, z, 'b'  )	% hits bottom only
                %                     case 3
                %                         plot3( x, y, z, 'g'  )	% hits surface only
                %                     otherwise
                %                         plot3( x, y, z, 'r' )
                %                 end
                
                set( gca, 'ZDir', 'Reverse' )   % plot with depth-axis positive down
                
                % update axis limits
                xmin = min( [ x xmin ] );
                xmax = max( [ x xmax xmin + .1 ] );

                ymin = min( [ y ymin ] );
                ymax = max( [ y ymax ymin + .1 ] );

                zmin = min( [ z zmin ] );
                zmax = max( [ z zmax zmin + .1 ] );
                if ( zmin == zmax ) % horizontal ray causes axis scaling problem
                    zmax = zmin + 1;
                end
                
                axis( [ xmin, xmax, ymin, ymax, zmin, zmax ] )
                
                if rem( ibeam, fix( Nalpha * Nbeta / 10 ) ) == 0    % flush graphics buffer every 10th ray
                    drawnow
                end
            end % next beam beta
        end	% next beam alpha
        
    end   % next Sy
end   % next Sx

fclose( fid );

drawnow
hold off
zoom on

if ( nargout == 1 )
    varargout{ 1 } = findobj( 'Type', 'Line' );   % return a handle to the lines in the figure
end

% fixed size for publications
if ( jkpsflag )
    set( gca, 'Units', 'centimeters' )
    set( gca, 'Position', [ 2 2 14.0  7.0 ] )
    set( gcf, 'PaperPositionMode', 'auto');
    
    %set( gcf, 'Units', 'centimeters' )
    %set( gcf, 'PaperPosition', [ 3 3 19.0 10.0 ] )
    set( gcf, 'Units', 'centimeters' )
    set( gcf, 'Position', [ 3 15 19.0 10.0 ] )
end