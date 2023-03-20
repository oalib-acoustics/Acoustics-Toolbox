function plotssp( envfil )
% plotssp.m
% Plots the sound speed profile
% useage:
%    plotssp( envfil )
%
% The envfil should be given without any extension
%
% MBP 03/2009

% read in the environmental file

[ ~, ~, ext ] = fileparts( envfil );

if ( ~strcmp( envfil, 'ENVFIL' ) && ~strcmp( ext, '.env' ) )
    envfil = [ envfil '.env' ]; % append extension
end

%[ ~, ~, SSP, Bdry, ~, ~, ~, ~, fid ] = read_env( envfil, 'BELLHOP' );
[ ~, ~, SSP, Bdry, fid ] = read_env_core( envfil );    % read in the environmental file

SSPType = Bdry.Top.Opt( 1 : 1 );

fclose( fid );

%figure
hold on

for medium = 1 : SSP.NMedia
   
   npts = round( SSP.raw( medium ).z(end) - SSP.raw( medium ).z( 1 ) ) + 1;
   if ( contains( SSPType, 'S' ) || contains( SSPType, 'P' ) )
      npts = 2 * npts;   % need to see internal points for cubic interpolation
   end
   
   z_eval = linspace( SSP.raw( medium ).z(1), SSP.raw( medium ).z( end ), npts );
   
   % plot the compression wave speed
   
   switch ( SSPType )
      case ( 'N' )   % n^2 Linear
         n2_eval = interp1( SSP.raw( medium ).z, 1.0 ./ ( SSP.raw( medium ).alphaR ).^2, z_eval, 'linear' );
         c_eval = real( 1.0 ./ sqrt( n2_eval ) );
      case ( 'P' )   % monotone PCHIP ACS (almost a cubic spline)
         c_eval = pchip_acs(  SSP.raw( medium ).z, real( SSP.raw( medium ).alphaR ), z_eval );
      case ( 'S' )   % Cubic Spline using not-a-knot boundary condition
         c_eval = spline( SSP.raw( medium ).z, real( SSP.raw( medium ).alphaR ), z_eval );
      otherwise      % piecewise Linear
         c_eval = real( interp1( SSP.raw( medium ).z, SSP.raw( medium ).alphaR, z_eval, 'linear' ) );
   end
   
   plot( real( SSP.raw( medium ).alphaR ), SSP.raw( medium ).z, 'ko' );
   hh = plot( c_eval, z_eval, 'b-' );
   set( hh, 'LineWidth', 2 );
   
   % plot the shear wave speed (if any non-zero values were supplied)
   
   if ( any( SSP.raw( medium ).betaR ) )
      
      switch ( SSPType )
         case ( 'N' )   % n^2 Linear
            n2_eval = interp1( SSP.raw( medium ).z, 1.0 ./ ( SSP.raw( medium ).betaR ).^2, z_eval, 'linear' );
            c_eval  = real( 1.0 ./ sqrt( n2_eval ) );
         case ( 'P' )   % monotone PCHIP ACS (almost a cubic spline)
            c_eval = pchip_acs(  SSP.raw( medium ).z, real( SSP.raw( medium ).betaR ), z_eval );
         case ( 'S' )   % Cubic Spline using not-a-knot boundary condition
            c_eval = spline( SSP.raw( medium ).z, real( SSP.raw( medium ).betaR ), z_eval );
         otherwise      % piecewise Linear
      end
      
      plot( real( SSP.raw( medium ).betaR ), SSP.raw( medium ).z, 'ko' );
      hh = plot( c_eval, z_eval, 'r-' );
      set( hh, 'LineWidth', 2 );
   end
   
end

set( gca, 'YDir', 'Reverse' )   % because view messes up the zoom feature
%axis IJ

xlabel( 'Sound Speed (m/s)' )
ylabel( 'Depth (m)' )
