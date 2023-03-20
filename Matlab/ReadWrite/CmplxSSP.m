function  [ SSP, Bdry ] = ComplexSSP( SSP, freq, AttenUnit )

% Set the imaginary part of the SSP based on the frequency
% note section at the end that has been commented out.
% This would need fixing for some models

%[ Bdry.Top.cp, Bdry.Top.cs, Bdry.Top.rho, ~, Bdry.Top.HS ] = topbot( fid, freq, Bdry.Top.BC, AttenUnit );

SSP.c   = [];
SSP.cs  = [];

% step through each point in the SSP
for medium = 1 : SSP.NMedia
   
   for ii = 1 : SSP.N( medium )
            
      cp = crci( SSP.raw( medium).alphaR( ii ), SSP.raw( medium).alphaI( ii ), freq, AttenUnit );
      cs = crci( SSP.raw( medium).betaR(  ii ), SSP.raw( medium).betaI(  ii ), freq, AttenUnit );
      SSP.c  = [ SSP.c;   cp   ];
      SSP.cs = [ SSP.cs;  cs   ];

   end
   
   if ( SSP.N( medium ) == 0 ) % calculate mesh automatically
      % choose a reference sound speed
      C = alphaR;
      if ( betaR > 0.0 ); C = betaR; end % shear?
      
      deltaz = 0.05 * C / freq;     % default sampling: 20 points per wavelength
      
      SSP.N( medium ) = round( ( SSP.depth( medium + 1 ) - SSP.depth( medium ) ) / deltaz );
      SSP.N( medium ) = max( SSP.N( medium ), 10 );     % require a minimum of 10 points
   end
   
end

% lower halfspace

Bdry.Bot.BC    = Bdry.Bot.Opt( 1: 1 );
[ Bdry.Bot.cp, Bdry.Bot.cs, Bdry.Bot.rho, ~, Bdry.Bot.HS ] = topbot( fid, freq, Bdry.Bot.BC, AttenUnit );

% Get rho, c just INSide the boundary (used later for reflection
% coefficients)
% I = NFirstAcoustic;
% Bdry.Top.rhoIns = SSP.rho( I );
% Bdry.Top.cIns   = SSP.c(   I );
% 
% I = Loc( NLastAcoustic ) + SSP.Npts( NLastAcoustic );
% Bdry.Bot.rhoIns = SSP.rho( I );
% Bdry.Bot.cIns   = SSP.c(   I );
