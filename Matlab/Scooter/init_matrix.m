function [ B1, B2, B3, B4, rho ] = init_matrix( NptsAll )

% Initializes arrays defining difference equations

global omega SSP Bdry
global NFirstAcoustic NLastAcoustic Loc H

SSPType = Bdry.Top.Opt( 1:1 );

% pre-allocate space for efficiency
B1  = zeros( NptsAll, 1 );
B2  = zeros( NptsAll, 1 );
B3  = zeros( NptsAll, 1 );
B4  = zeros( NptsAll, 1 );
rho = zeros( NptsAll, 1 );

cMin           = 1.0E6;
NFirstAcoustic = 0;
Loc( 1 )       = 0;
omega2 = omega^2;

% *** Loop over media ***

for Medium = 1:SSP.NMedia
   if ( Medium ~= 1 )
      Loc( Medium ) = Loc( Medium - 1 ) + SSP.N( Medium - 1 ) + 1;
   end
   N1   = SSP.N(   Medium ) + 1;
   ii   = Loc( Medium ) + 1;
   
   % calculate indices in user specified SSP for that medium
   if Medium == 1
      I1 = 1;
   else
      I1 = 1 + sum( SSP.Npts( 1:Medium-1) );
   end
   I2 = I1 + SSP.Npts( Medium ) - 1;
   II = I1 : I2;
   ZT  = linspace( SSP.z( I1 ), SSP.z( I2 ), N1 );
   jj = ii:ii + SSP.N( Medium );
   
   switch ( SSPType )
      case ( 'N' )   % n2-linear approximation to SSP
         cp        = 1./ sqrt( interp1( SSP.z( II ), 1./SSP.c(  II ).^2,  ZT ) );
         if ( any( SSP.cs( II ) ) )
            cs     = 1./ sqrt( interp1( SSP.z( II ), 1./SSP.cs( II ).^2,  ZT ) );
         else
            cs = zeros( size( SSP.cs ) );
         end
         rho( jj ) = interp1( SSP.z( II ), SSP.rho(   II ),  ZT );
      case ( 'C' )   % C-LINEAR approximation to SSP
         cp        = interp1( SSP.z( II ), SSP.c(   II ),  ZT );
         cs        = interp1( SSP.z( II ), SSP.cs(  II ),  ZT );
         rho( jj ) = interp1( SSP.z( II ), SSP.rho( II ),  ZT );
      case ( 'P' )   % monotone PCHIP ACS (almost a cubic spline)
         cp        = pchip_acs( SSP.z( II ), SSP.c(   II ),  ZT );
         cs        = pchip_acs( SSP.z( II ), SSP.cs(  II ),  ZT );
         rho( jj ) = pchip_acs( SSP.z( II ), SSP.rho( II ),  ZT );
      case ( 'S' )   % Cubic Spline using not-a-knot boundary condition
         cp        = spline( SSP.z( II ), SSP.c(   II ),  ZT );
         cs        = spline( SSP.z( II ), SSP.cs(  II ),  ZT );
         rho( jj ) = spline( SSP.z( II ), SSP.rho( II ),  ZT );
      case ( 'A' )
         error( '    ANALYTIC SSP option not available' )
   end

   % note: calculations here for first, last acoustic media redundant (see
   % readenv)
   if ( ~any( cs ) )   % *** Case of an acoustic Medium ***  
      %Mater( Medium, : )  =  'ACOUSTIC';
      if ( NFirstAcoustic == 0 )
         NFirstAcoustic = Medium;
      end
      NLastAcoustic  = Medium;
      
      cMinV = min( real( cp ) );
      cMin  = min( cMin, cMinV );
      B1( jj ) = omega2 ./ cp.^2;
      
   else                % *** Case of an elastic medium ***
      %Mater( Medium, : ) = 'ELASTIC ';
      TwoH = 2.0 * H( Medium );
      
         cMinV = min( real( cs ) );
         cMin  = min( cMin, cMinV );

         cp2 = cp.^2;
         cs2 = cs.^2;

         B1(  jj ) = TwoH ./ ( rho( jj )' .* cs2 );
         B2(  jj ) = TwoH ./ ( rho( jj )' .* cp2 );
         B3(  jj ) = 4.0 * TwoH * rho( jj )' .* cs2 .* ( cp2 - cs2 ) ./ cp2;
         B4(  jj ) = TwoH * ( cp2 - 2.0 .* cs2 ) ./ cp2;
         rho( jj ) = TwoH * omega2 * rho( jj )';
   end
end % next Medium
