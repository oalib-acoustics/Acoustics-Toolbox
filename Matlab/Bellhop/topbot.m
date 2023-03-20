function [ cp, cs, rho, HS ] = topbot( fid, freq, BCType, AttenUnit )

% Handles top and bottom boundary conditions

% Input:
%   ENVFIL: Environmental file
%   freq:   frequency
%   BCType: Boundary condition type

% Output:
%  Bdry.cp:    P-wave speed in halfspace
%  Bdry.cs:    S-wave speed in halfspace
%  Bdry.rho:   density in halfspace

global alphaR betaR rhoR alphaI betaI

% *** Echo to PRTFIL user's choice of boundary condition ***

switch ( BCType )
case ( 'V' )
   disp( '    VACUUM' )
case ( 'R' )
   disp( '    Perfectly RIGID' )
case ( 'A' )
   disp( '    ACOUSTO-ELASTIC half-space' )
case ( 'F' )
   disp( '    FILE used for reflection loss' )
case ( 'W' ) 
   disp( '    Writing an IRC file' )
case ( 'P' )
   disp( '    reading PRECALCULATED IRC' )
otherwise
   error( 'Fatal error: Unknown boundary condition type' )
end

% ****** Read in BC parameters depending on particular choice ******

cp  = 0.0;
cs  = 0.0;
rho = 0.0;

HS = [];

% *** Half-space properties ***

if ( BCType == 'A' )
   ztmp   = fscanf( fid, '%f', 1 );
   alphaRtemp = fscanf( fid, '%f', 1 );
   betaRtemp  = fscanf( fid, '%f', 1 );
   rhoRtemp   = fscanf( fid, '%f', 1 );
   alphaItemp = fscanf( fid, '%f', 1 );
   betaItemp  = fscanf( fid, '%f', 1 );
   fgetl( fid );
   % if values read in copy over, otherwise use defaults
   if (~isempty( alphaRtemp ) ); alphaR = alphaRtemp; end
   if (~isempty( betaRtemp  ) ); betaR  = betaRtemp;  end
   if (~isempty( rhoRtemp   ) ); rhoR   = rhoRtemp;   end
   if (~isempty( alphaItemp ) ); alphaI = alphaItemp; end
   if (~isempty( betaItemp  ) ); betaI  = betaItemp;  end

   fprintf( '%10.2f    %10.2f    %10.2f    %10.2f    %10.4f    %10.4f \n', ztmp, alphaR, betaR, rhoR, alphaI, betaI )

   cp  = crci( alphaR, alphaI, freq, AttenUnit );
   cs  = crci( betaR,  betaI,  freq, AttenUnit );
   rho = rhoR;
   
   % raw values in user units
   HS.alphaR = alphaR;
   HS.alphaI = alphaI;
   HS.betaR  = betaR;
   HS.betaI  = betaI;
   HS.rho    = rhoR;
end
