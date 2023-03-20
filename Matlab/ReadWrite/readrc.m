function readrc( trcfil, brcfil, TopRC, BotRC )

% Read the top and bottom reflection coefficients

global thetaBot RBot phiBot NBotPts thetaTop RTop phiTop NTopPts


%%
if ( TopRC == 'F' )
    
    disp( 'Reading top    reflection coefficient' )
    
    if ( strcmp( trcfil, 'TRCFIL' ) == 0 &&  ~contains( trcfil, '.trc' )  )
        trcfil = [ trcfil '.trc' ]; % append extension
    end
    
    fid = fopen( trcfil, 'r' );
    if ( fid == -1 )
        warndlg( 'Top Refl. Coef. File, TRCFIL. does not exist', 'Warning' )
    end
    
    NTopPts  = fscanf( fid, '%i', 1 );
    fprintf( 'Number of points in top    reflection coefficent = %i \n', NTopPts )
    thetaTop = zeros( NTopPts, 1 );
    RTop     = zeros( NTopPts, 1 );
    phiTop   = zeros( NTopPts, 1 );

    for ii = 1 : NTopPts
        thetaTop( ii ) = fscanf( fid, '%f', 1 );
        RTop(     ii ) = fscanf( fid, '%f', 1 );
        phiTop(   ii ) = fscanf( fid, '%f', 1 );
    end
    fclose( fid );
    phiTop = pi / 180 * phiTop;   % convert to radians
    
else
    thetaTop = 0;
    RTop     = 0;
    phiTop   = 0;
    NTopPts  = 0;
end

validateattributes( thetaTop, { 'numeric' }, { 'vector', 'nondecreasing' } )

%%

if ( BotRC == 'F' )
    disp( 'Reading bottom reflection coefficient' )
    
    if ( strcmp( brcfil, 'BRCFIL' ) == 0 &&  ~contains( brcfil, '.brc' )  )
        brcfil = [ brcfil '.brc' ]; % append extension
    end
    
    fid = fopen( brcfil, 'r' );
    if ( fid == -1 )
        warndlg( 'Bot Refl. Coef. File, BRCFIL. does not exist', 'Warning' )
    end
    
    NBotPts  = fscanf( fid, '%i', 1 );
    fprintf( 'Number of points in bottom reflection coefficent = %i \n', NBotPts )
    thetaBot = zeros( NBotPts, 1 );
    RBot     = zeros( NBotPts, 1 );
    phiBot   = zeros( NBotPts, 1 );

    for ii = 1 : NBotPts
        thetaBot( ii ) = fscanf( fid, '%f', 1 );
        RBot(     ii ) = fscanf( fid, '%f', 1 );
        phiBot(   ii ) = fscanf( fid, '%f', 1 );
    end
    fclose( fid );
    phiBot = pi / 180 * phiBot;   % convert to radians
    
else
    thetaBot = 0;
    RBot     = 0;
    phiBot   = 0;
    NBotPts  = 0;
end

validateattributes( thetaBot, { 'numeric' }, { 'vector', 'nondecreasing' } )

