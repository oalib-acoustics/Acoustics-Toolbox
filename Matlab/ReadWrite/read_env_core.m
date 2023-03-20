function  [ TitleEnv, freq, SSP, Bdry, fid ] = read_env_core( envfil )
% Read the core of the environmental file
% This is the part used by all the models in the Acoustics Toolbox

global HV NFirstAcoustic NLastAcoustic
global alphaR betaR rhoR alphaI betaI
global T S pH z_bar

alphaR = 1500;   % defaults
betaR  = 0;
rhoR   = 1;
alphaI = 0;
betaI  = 0;

NFirstAcoustic = 0;

Bdry.Top.cp    = 0.0;
Bdry.Top.cs    = 0.0;
Bdry.Top.rho   = 0.0;
Bdry.Bot.cp    = 2000.0;
Bdry.Bot.cs    = 0.0;
Bdry.Bot.rho   = 2.0;

fid = fopen( envfil, 'r' );
if ( fid == -1 )
   error( 'Unable to open environmental file', 'readenv' );
end

TitleEnv = fgetl( fid );
% Extract letters between the quotes
nchars   = strfind( TitleEnv, '''' );   % find quotes
TitleEnv = TitleEnv( nchars( 1 ) + 1 : nchars( 2 ) - 1 );

disp( TitleEnv )

freq     = fscanf( fid, '%f', 1 );
fprintf( 'Frequency = %d Hz \n', freq )
fgetl( fid );

SSP.NMedia   = fscanf( fid, '%i', 1 );
fprintf( 'Number of media = %i \n\n', SSP.NMedia )
fgetl( fid );

TopOpt   = fgetl( fid );

% Extract option letters between the quotes
nchars = strfind( TopOpt, '''' );   % find quotes
Bdry.Top.Opt   = [ TopOpt( nchars( 1 ) + 1 : nchars( 2 ) - 1 ) blanks( 7 - ( nchars( 2 ) - nchars( 1 ) ) ) ];

% convert the deprecated '*' option to '~'
if ( Bdry.Top.Opt( 5 : 5 ) == '*' )
   Bdry.Top.Opt( 5 : 5 ) = '~';
end

SSPType(   1 : 1 ) = Bdry.Top.Opt( 1 : 1 );
Bdry.Top.BC        = Bdry.Top.Opt( 2 : 2 );
AttenUnit( 1 : 2 ) = Bdry.Top.Opt( 3 : 4 );

%     *** SSP approximation options ***

switch ( SSPType )
   case ( 'N' )
      disp( '    N2-Linear approximation to SSP' )
   case ( 'C' )
      disp( '    C-Linear approximation to SSP' )
   case ( 'P' )
      disp( '    PCHIP approximation to SSP' )
   case ( 'S' )
      disp( '    Spline approximation to SSP' )
   case ( 'Q' )
      disp( '    Quadrilateral approximation to range-dependent SSP' )
   case ( 'H' )
      disp( '    Hexahedral approximation to range and depth dependent SSP' )
   case ( 'A' )
      disp( '    Analytic SSP option' )
   otherwise
      fclose all;
      error( 'Fatal error: Unknown option for SSP approximation' )
end

%     *** Attenuation options ***

switch ( AttenUnit( 1 : 1 ) )
   case ( 'N' )
      disp( '    Attenuation units: nepers/m' )
   case ( 'F' )
      disp( '    Attenuation units: dB/mkHz' )
   case ( 'M' )
      disp( '    Attenuation units: dB/m' )
   case ( 'W' )
      disp( '    Attenuation units: dB/wavelength' )
   case ( 'Q' )
      disp( '    Attenuation units: Q' )
   case ( 'L' )
      disp( '    Attenuation units: Loss tangent' )
   otherwise
      fclose all;
      error( 'Fatal error: Unknown attenuation units' )
end

%     *** optional addition of volume attenuation using standard formulas

if ( length( Bdry.Top.Opt ) >= 4 )
   switch ( Bdry.Top.Opt(4:4) )
      case ( 'T' )
         disp( '    THORP attenuation added' )
      case ( 'F' )
         disp( '    Francois-Garrison attenuation added' )
         SSP.T     = fscanf( fid, '%f', 1 );
         SSP.S     = fscanf( fid, '%f', 1 );
         SSP.pH    = fscanf( fid, '%f', 1 );
         SSP.z_bar = fscanf( fid, '%f', 1 );
         fgetl( fid );
         
         T     = SSP.T;
         S     = SSP.S;
         pH    = SSP.pH;
         z_bar = SSP.z_bar;
         
         fprintf(  '        T =  %4.1f degrees   S = %4.1f psu   pH = %4.1f   z_bar = %6.1f m \n', ...
            SSP.T, SSP.S, SSP.pH, SSP.z_bar );
   end
end

if ( length( Bdry.Top.Opt ) >= 5 )
   switch ( Bdry.Top.Opt(5:5) )
      case ( '*' )
         disp( '    Development options enabled' )
   end
end

[ Bdry.Top.cp, Bdry.Top.cs, Bdry.Top.rho, Bdry.Top.HS ] = topbot( fid, freq, Bdry.Top.BC, AttenUnit );

%% main loop to readin in SSP
fprintf( '\n       z          alphaR         betaR           rho        alphaI         betaI' );
fprintf( '\n      (m)          (m/s)         (m/s)         (g/cm^3)      (m/s)         (m/s) \n' );
SSP.z   = [];
SSP.c   = [];
SSP.cs  = [];
SSP.rho = [];
for medium = 1 : SSP.NMedia
   if ( medium == 1 )
      Loc( medium ) = 0;
   else
      Loc( medium ) = Loc( medium - 1 ) + SSP.Npts( medium - 1 );
   end
   
   % number of points in finite-difference grid
   SSP.N(     medium  ) = fscanf( fid, '%i', 1 );
   SSP.sigma( medium  ) = fscanf( fid, '%f', 1 );
   SSP.depth( medium+1) = fscanf( fid, '%f', 1 );
   
   
   fprintf( '    ( Number of points = %i  Roughness = %6.2f  Depth = %8.2f ) \n', ...
      SSP.N( medium ), SSP.sigma( medium ), SSP.depth( medium+1 ) );
   fgetl( fid );
   
   % read in the SSP
   
   for ii = 1 : 9999999   % inf generated an error in latest Matlab
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
      
      cp = crci( alphaR, alphaI, freq, AttenUnit );
      cs = crci( betaR,  betaI,  freq, AttenUnit );
      SSP.z   = [ SSP.z;   ztmp ];   % add in to existing vector
      SSP.c   = [ SSP.c;   cp   ];
      SSP.cs  = [ SSP.cs;  cs   ];
      SSP.rho = [ SSP.rho; rhoR ];
      
      SSP.raw( medium).z( ii )      = ztmp;   % add in to existing vector
      SSP.raw( medium).alphaR( ii ) = alphaR;
      SSP.raw( medium).alphaI( ii ) = alphaI;
      SSP.raw( medium).betaR( ii )  = betaR;
      SSP.raw( medium).betaI( ii )  = betaI;
      SSP.raw( medium).rho( ii )    = rhoR;
      
      % check for end of this layer
      if ( ztmp == SSP.depth( medium+1 ) )
         SSP.depth( 1 ) = SSP.z( 1 );
         break
      end
   end
   
   if ( SSP.N( medium ) == 0 ) % calculate mesh automatically
      % choose a reference sound speed
      C = alphaR;
      if ( betaR > 0.0 ); C = betaR; end % shear?
      
      deltaz = 0.05 * C / freq;     % default sampling: 20 points per wavelength
      
      SSP.N( medium ) = round( ( SSP.depth( medium + 1 ) - SSP.depth( medium ) ) / deltaz );
      SSP.N( medium ) = max( SSP.N( medium ), 10 );     % require a minimum of 10 points
   end
   
   fprintf( '    Number of points = %i \n', SSP.N( medium ) );
   
   % keep track of first and last acoustic medium
   if ( ~any( cs ) )   % shear anywhere?
      if ( NFirstAcoustic == 0 )
         NFirstAcoustic = medium;
      end
      NLastAcoustic  = medium;
   end
   
   % stuff for Bellhop
   if ( medium == 1 )
      HV       = diff( SSP.z ); % layer thicknesses
      SSP.cz   = diff( SSP.c ) ./ HV; % gradient of ssp (centered-difference approximation)
   end
   SSP.Npts( medium ) = ii;
   if ( medium == 1 ); SSP.depth( 1 ) = SSP.z( 1 ); end
end

%% lower halfspace

BotOpt = fgetl( fid );
nchars = strfind( BotOpt, '''' );   % find quotes

% Extract option letters between the quotes
Bdry.Bot.Opt   = [ BotOpt( nchars( 1 ) + 1 : nchars( 2 ) - 1 ) blanks( 3 - ( nchars( 2 ) - nchars( 1 ) ) ) ];

% convert the deprecated '*' option to '~'
if ( Bdry.Bot.Opt( 2 : 2 ) == '*' )
   Bdry.Bot.Opt( 2 : 2 ) = '~';
end

Bdry.Bot.BC    = Bdry.Bot.Opt( 1: 1 );
[ Bdry.Bot.cp, Bdry.Bot.cs, Bdry.Bot.rho, Bdry.Bot.HS ] = topbot( fid, freq, Bdry.Bot.BC, AttenUnit );

Bdry.Top.depth = SSP.depth( 1 );
Bdry.Bot.depth = SSP.depth( SSP.NMedia + 1 );

% Get rho, c just INSide the boundary (used later for reflection
% coefficients)
I = NFirstAcoustic;
Bdry.Top.rhoIns = SSP.rho( I );
Bdry.Top.cIns   = SSP.c(   I );

I = Loc( NLastAcoustic ) + SSP.Npts( NLastAcoustic );
Bdry.Bot.rhoIns = SSP.rho( I );
Bdry.Bot.cIns   = SSP.c(   I );
