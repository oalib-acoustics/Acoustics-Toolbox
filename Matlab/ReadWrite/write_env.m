function write_env( envfil, model, TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, varargin )

% Write an environmental file
% mbp 2009
% note: I had an unusual case with a parabolic mirror where round-off in
% the receiver ranges was a problem (+/- 1 m in range)
% Changed the output format from %6.2f to %6f to accommodate that.
% Similar issues may occur elsewhere in the code below ...

if ( strcmp( envfil, 'ENVFIL' ) == 0 && ~strcmp( envfil( end-3: end ), '.env' ) )
  envfil = [ envfil '.env' ]; % append extension
end

if ( size( varargin ) == 0 )
    fid = fopen( envfil, 'wt' );   % create new envfil
else
    fid = fopen( envfil, 'at' );   % append to existing envfil
end

if ( fid == -1 )
    disp( envfil )
    error( 'Unable to create environmental file', 'write_env' );
end

model = upper( model );   % convert to uppercase

fprintf( fid, '''%s'' ! Title \n', TitleEnv );
fprintf( fid, '%8.2f  \t \t \t ! Frequency (Hz) \n', freq );
fprintf( fid, '%5i    \t \t \t ! NMedia \n', SSP.NMedia );
fprintf( fid, '''%s'' \t \t \t ! Top Option \n', Bdry.Top.Opt );

if ( Bdry.Top.Opt( 2:2 ) == 'A' )
    fprintf( fid, '    %6.2f %6.2f %6.2f %6.2g %6.2f %6.2f /  \t ! upper halfspace \n', SSP.depth( 1 ), ...
        Bdry.Top.HS.alphaR, Bdry.Top.HS.betaR, Bdry.Top.HS.rho, Bdry.Top.HS.alphaI, Bdry.Top.HS.betaI );
end

% SSP
for medium = 1 : SSP.NMedia
    
    fprintf( fid, '%5i %4.2f %6.2f \t ! N sigma depth \n', SSP.N( medium ), SSP.sigma( medium ), SSP.depth( medium+1 ) );
    for ii = 1 : length( SSP.raw( medium ).z )
        fprintf( fid, '\t %6.2f %6.2f %6.2f %6.2g %10.6f %6.2f / \t ! z c cs rho \n', ...
            [ SSP.raw( medium ).z( ii ) ...
              SSP.raw( medium ).alphaR( ii ) SSP.raw( medium ).betaR( ii ) SSP.raw( medium ).rho( ii ) ...
              SSP.raw( medium ).alphaI( ii ) SSP.raw( medium ).betaI( ii ) ] );
    end
end

% lower halfspace
fprintf( fid, '''%s'' %6.2f  \t \t ! Bottom Option, sigma \n', Bdry.Bot.Opt, 0.0 ); % SSP.sigma( 2 ) );

if ( Bdry.Bot.Opt( 1:1 ) == 'A' )
    fprintf( fid, '    %6.2f %6.2f %6.2f %6.2g %6.2f %6.2f /  \t ! lower halfspace \n', SSP.depth( SSP.NMedia+1 ), ...
        Bdry.Bot.HS.alphaR, Bdry.Bot.HS.betaR, Bdry.Bot.HS.rho, Bdry.Bot.HS.alphaI, Bdry.Bot.HS.betaI );
end

if( strmatch( model, strvcat( 'SCOOTER', 'KRAKEN', 'KRAKENC', 'SPARC' ), 'exact' ) )
    fprintf( fid, '%6.0f %6.0f \t \t ! cLow cHigh (m/s) \n', cInt.Low, cInt.High );   % phase speed limits
    fprintf( fid, '%8.2f \t \t \t ! RMax (km) \n', RMax );    % maximum range
end

% source depths

fprintf( fid, '%5i \t \t \t \t ! NSz \n', length( Pos.s.z ) );

if ( length( Pos.s.z ) >= 2 && equally_spaced( Pos.s.z ) )
    fprintf( fid, '    %6f %6f', Pos.s.z( 1 ), Pos.s.z( end ) );
else
    fprintf( fid, '    %6f  ', Pos.s.z );
end

fprintf( fid, '/ \t ! Sz(1)  ... (m) \n' );

% receiver depths

fprintf( fid, '%5i \t \t \t \t ! NRz \n', length( Pos.r.z ) );

if ( length( Pos.r.z ) >= 2 && equally_spaced( Pos.r.z ) )
    fprintf( fid, '    %6f %6f ', Pos.r.z( 1 ), Pos.r.z( end ) );
else
    fprintf( fid, '    %6f  ', Pos.r.z );
end

fprintf( fid, '/ \t ! Rz(1)  ... (m) \n' );

% receiver ranges
if ( strcmp( model, 'BELLHOP' ) ||  strcmp( model, 'simplePE' ) )
    fprintf( fid, '%5i \t \t \t \t ! NRr \n', length( Pos.r.r ) );
    
    if ( length( Pos.r.r ) >= 2 && equally_spaced( Pos.r.r ) )
        fprintf( fid, '    %6f %6f', Pos.r.r( 1 ), Pos.r.r( end ) );
    else
        fprintf( fid, '    %6f ', Pos.r.r );
    end
    fprintf( fid, '/ \t ! Rr(1)  ... (km) \n' );
    write_bell( fid, Beam );
end

fclose( fid );
