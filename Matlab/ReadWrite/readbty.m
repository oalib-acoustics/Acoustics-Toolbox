function readbty( btyfil, BotBTY, depthB, rBox )

%READBTY Read a  btyfil into the workspace variables
% xBot contains the range depth info with ranges converted to meters


global Bdry
global xBot tBot nBot tBotNode nBotNode RLenBot kappaBot NbtyPts btyType

if ( BotBTY == '*' || BotBTY == '~' )
   fprintf( '\n_______________________ \n' )
   disp( 'Using bottom-bathymetry file' )
   
   if ( strcmp( btyfil, 'BTYFIL' ) == 0 &&  ~contains( btyfil, '.bty' )  )
      btyfil = [ btyfil '.bty' ]; % append extension
   end
   
   fid = fopen( btyfil, 'r' );
   if ( fid == -1 )
      error( 'Bathymetry file does not exist' )
   end
   
   btyType = fgetl( fid );
   
   % Extract option letter between the quotes
   nchars = strfind( btyType, '''' );   % find quotes
   btyType = [ btyType( nchars( 1 ) + 1 : nchars( 2 ) - 1 ) blanks( 2 - ( nchars( 2 ) - nchars( 1 ) ) ) ];
   
   switch ( btyType( 1 : 1 ) )
      case ( 'L' )
         disp( 'Piecewise-linear approximation to bathymetry' )
      case ( 'C' )
         disp( 'Curvilinear approximation to bathymetry' )
      otherwise
         fclose all;
         disp( btyType )
         error( 'Fatal error: Unknown option for bathymetry type' )
   end
   
   NbtyPts = fscanf( fid, '%i', 1 );
   fprintf( 'Number of bathymetry points = %i \n\n', NbtyPts )
   fprintf( ' Range (km)     Depth (m) \n' )
   
   % btymat  = fscanf( fid, '%f', [ 2 NbtyPts ] );
   
   xBot = zeros( 2, NbtyPts + 2 );
   
   for ii = 1 : NbtyPts
      xBot( :, ii+1 ) = fscanf( fid, '%f', 2 );
      if ( ii == NbtyPts && ii > 11 )
         disp( '    ...' )
      end
      
      if ( ii < 11 || ii == NbtyPts )   % echo up to 11 values
         fprintf( '%9.5g    %9.5g \n', xBot( :, ii+1 ) );
      end
      
      if ( xBot( 2, ii ) > depthB )
         error( 'Bathymetry goes below last tabulated SSP value' )
      end
   end
   fclose( fid );  % close the bathymetry file

   xBot( 1, : ) = 1000.0 * xBot( 1, : );   % Convert ranges in km to m
   
   % extend the bathymetry to +/- infinity in a piecewise constant fashion
   NbtyPts = NbtyPts + 2;
   xBot( 1, 1 )       = -1e50;
   xBot( 2, 1 )       = xBot( 2, 2);
   xBot( 1, NbtyPts ) = +1e50;
   xBot( 2, NbtyPts ) = xBot( 2, NbtyPts - 1 );
   
   % compute tangent and outward-pointing normal to Bot
   tBot = zeros( 2, NbtyPts-1 );   % in case tBot exists from a previous, different run
   nBot = zeros( 2, NbtyPts-1 );
   
   tBot( 1, : ) = xBot( 1, 2:NbtyPts ) - xBot( 1, 1:NbtyPts - 1 );
   tBot( 2, : ) = xBot( 2, 2:NbtyPts ) - xBot( 2, 1:NbtyPts - 1 );
   
   RLenBot = sqrt( tBot( 1, : ).^ 2 + tBot( 2, : ).^ 2 );
   
   tBot( 1, : ) = tBot( 1, : ) ./ RLenBot;
   tBot( 2, : ) = tBot( 2, : ) ./ RLenBot;
   
   nBot( 1, : ) = -tBot( 2, : );
   nBot( 2, : ) =  tBot( 1, : );

   if ( btyType( 1 : 1 ) == 'C' ) % Curvilinear option: compute tangent and normal at node by averaging normals on adjacent segments
      
      tBotNode = zeros( 2, NbtyPts );
      nBotNode = zeros( 2, NbtyPts );
      for ii = 2 : NbtyPts - 1
         tBotNode( :, ii ) = 0.5 * ( tBot( :, ii - 1 ) + tBot( :, ii ) );
         nBotNode( :, ii ) = 0.5 * ( nBot( :, ii - 1 ) + nBot( :, ii ) );
      end
      
      % end points
      tBotNode( :, 1       ) = [ 1.0, 0.0 ];   % tangent left-end  node
      tBotNode( :, NbtyPts ) = [ 1.0, 0.0 ];   % tangent right-end node
      nBotNode( 1, : ) = -tBotNode( 2, : );
      nBotNode( 2, : ) = +tBotNode( 1, : );
      
      % compute curvature in each segment
      
      phi = atan2( tBotNode( 2, : ), tBotNode( 1, : ) )';
      kappaBot = diff( phi( 1 : end ) ) ./ RLenBot( 1 : end )'; % this is curvature = dphi/ds
   else
      kappaBot = zeros( NbtyPts + 1, 1 );
   end
   
else   % no bathymetry given, use SSP depth for flat Bot
   NbtyPts = 3;
   xBot = zeros( 2, NbtyPts );
   
   xBot( :, 1 ) = [ -1000 * rBox; Bdry.Bot.depth ];
   xBot( :, 2 ) = [            0; Bdry.Bot.depth ];
   xBot( :, 3 ) = [  1000 * rBox; Bdry.Bot.depth ];
   
   tBot = zeros( 2, NbtyPts - 1 );   % in case tBot exists from a previous, different run
   nBot = zeros( 2, NbtyPts - 1 );
   
   tBot( :, 1 ) = [ 1.0;  0.0 ];   % tangent to Bot
   nBot( :, 1 ) = [ 0.0;  1.0 ];   % outward-pointing normal
   tBot( :, 2 ) = [ 1.0;  0.0 ];   % tangent to Bot
   nBot( :, 2 ) = [ 0.0;  1.0 ];   % outward-pointing normal
   tBot( :, 3 ) = [ 1.0;  0.0 ];   % tangent to Bot
   nBot( :, 3 ) = [ 0.0;  1.0 ];   % outward-pointing normal
   
   RLenBot = sqrt( tBot( 1, : ).^ 2 + tBot( 2, : ).^ 2 );
   
   kappaBot( 1, 1 ) = 0;
   kappaBot( 2, 1 ) = 0;
   
end
