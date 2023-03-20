function readati( atifil, TopATI, depthT, rBox )

% Read an atifil into the workspace variables

global Bdry
global xTop tTop nTop tTopNode nTopNode RLenTop kappaTop NatiPts atiType

if ( TopATI == '*' || TopATI == '~' )
   fprintf( '\n_______________________ \n' )
   disp( 'Using top-altimetry file' )
   
   if ( strcmp( atifil, 'ATIFIL' ) == 0 &&  ~contains( atifil, '.ati' )  )
      atifil = [ atifil '.ati' ]; % append extension
   end
   
   fid = fopen( atifil, 'r' );
   if ( fid == -1 )
      error( 'ATIFIL does not exist' )
   end
   
   atiType = fgetl( fid );
   
   % Extract option letter between the quotes
   nchars = strfind( atiType, '''' );   % find quotes
   atiType = [ atiType( nchars( 1 ) + 1 : nchars( 2 ) - 1 ) blanks( 2 - ( nchars( 2 ) - nchars( 1 ) ) ) ];
   
   switch ( atiType )
      case ( 'L' )
         disp( 'Piecewise-linear approximation to altimetry' )
      case ( 'C' )
         disp( 'Curvilinear approximation to altimetry' )
      otherwise
         fclose all;
         disp( atiType )
         error( 'Fatal error: Unknown option for altimetry type' )
   end
   
   % [ xTop, NatiPts ] = readvector( fid );
   NatiPts = fscanf( fid, '%i', 1 );
   fprintf( 'Number of altimetry points = %i \n\n', NatiPts )
   fprintf( ' Range (km)     Depth (m) \n' )
   
   %atimat  = fscanf( fid, '%f', [ 2 Natipts ] );
   
   xTop = zeros( 2, NatiPts+2 );
   
   for ii = 1 : NatiPts
      xTop( :, ii+1 ) = fscanf( fid, '%f', 2 );
      
      if ( ii == NatiPts && ii > 11 )
         disp( '    ...' )
      end
      
      if ( ii < 11 || ii == NatiPts )   % echo up to 11 values
         fprintf( '%9.5g    %9.5g \n', xTop( :, ii+1 ) );
      end
      if ( xTop( 2, ii ) < depthT )
         error( 'Altimetry goes above first tabulated SSP value' )
      end
   end
   fclose( fid );   % close the altimetry file
   
   xTop( 1, : ) = 1000.0 * xTop( 1, : );   % Convert ranges in km to m
   
   % extend the bathymetry to +/- infinity in a piecewise constant fashion
   NatiPts = NatiPts + 2;
   xTop( 1, 1 )       = -1e50;
   xTop( 2, 1 )       = xTop( 2, 2);
   xTop( 1, NatiPts ) = +1e50;
   xTop( 2, NatiPts ) = xTop( 2, NatiPts - 1 );
   
   % compute tangent and outward-pointing normal to top
   tTop = zeros( 2, NatiPts-1 );   % in case tTop exists from a previous, different run
   nTop = zeros( 2, NatiPts-1 );
   
   tTop( 1, : ) = xTop( 1, 2:NatiPts ) - xTop( 1, 1:NatiPts - 1 );
   tTop( 2, : ) = xTop( 2, 2:NatiPts ) - xTop( 2, 1:NatiPts - 1 );
   
   RLenTop = sqrt( tTop( 1, : ).^ 2 + tTop( 2, : ).^ 2 );
   
   tTop( 1, : ) = tTop( 1, : ) ./ RLenTop;
   tTop( 2, : ) = tTop( 2, : ) ./ RLenTop;
   
   nTop( 1, : ) =  tTop( 2, : );
   nTop( 2, : ) = -tTop( 1, : );
   
   if ( atiType( 1 : 1 ) == 'C' ) % Curvilinear option: compute tangent and normal at node by averaging normals on adjacent segments
      
      tTopNode = zeros( 2, NatiPts );
      nTopNode = zeros( 2, NatiPts );
      for ii = 2 : NatiPts - 1
         tTopNode( :, ii ) = 0.5 * ( tTop( :, ii - 1 ) + tTop( :, ii ) );
         nTopNode( :, ii ) = 0.5 * ( nTop( :, ii - 1 ) + nTop( :, ii ) );
      end
      
      % end points
      tTopNode( :, 1       ) = [ 1.0, 0.0 ];   % tangent left-end  node
      tTopNode( :, NatiPts ) = [ 1.0, 0.0 ];   % tangent right-end node
      nTopNode( 1, : ) = +tTopNode( 2, : );
      nTopNode( 2, : ) = -tTopNode( 1, : );
      
      % compute curvature in each segment
      
      phi = atan2( tTopNode( 2, : ), tTopNode( 1, : ) )';
      kappaTop = diff( phi( 1 : end ) ) ./ RLenTop( 1 : end )'; % this is curvature = dphi/ds
   else
      kappaTop = zeros( NatiPts + 1, 1 );
   end
else   % no bathymetry given, use SSP depth for flat top
   NatiPts = 3;
   xTop = zeros( 2, NatiPts );
   
   xTop( :, 1 ) = [ -1000 * rBox; Bdry.Top.depth ];
   xTop( :, 2 ) = [     0; Bdry.Top.depth ];
   xTop( :, 3 ) = [  1000 * rBox; Bdry.Top.depth ];
   
   tTop = zeros( 2, NatiPts-1 );   % in case tTop exists from a previous, different run
   nTop = zeros( 2, NatiPts-1 );
   
   tTop( :, 1 ) = [ 1.0;  0.0 ];   % tangent to top
   nTop( :, 1 ) = [ 0.0; -1.0 ];   % outward-pointing normal
   tTop( :, 2 ) = [ 1.0;  0.0 ];   % tangent to top
   nTop( :, 2 ) = [ 0.0; -1.0 ];   % outward-pointing normal
   tTop( :, 3 ) = [ 1.0;  0.0 ];   % tangent to top
   nTop( :, 3 ) = [ 0.0; -1.0 ];   % outward-pointing normal
   
   RLenTop = sqrt( tTop( 1, : ).^ 2 + tTop( 2, : ).^ 2 );
   
   kappaTop( 1, 1 ) = 0;
   kappaTop( 2, 1 ) = 0;
   
end