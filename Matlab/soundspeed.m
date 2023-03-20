function ssp=soundspeed(S,T,D,equation)
% SOUNDSPEED Various sound-speed equations.
%
%          1)  SOUNDSPEED(S,T,D) returns the sound speed (m/sec) given vectors
%             of salinity (ppt), temperature (deg C) and DEPTH (m) using
%             the formula of Mackenzie:
%
%           Mackenzie, K.V. "Nine-term Equation for Sound Speed in the Oceans",
%           J. Acoust. Soc. Am. 70 (1981), 807-812.
%
%
%          2) SOUNDSPEED(S,T,P,'del grosso') returns the sound speed (m/sec) 
%            given vectors of salinity (ppt), temperature (deg C), and 
%            PRESSURE (dbar) using the Del Grosso equation:
%
%            Del Grosso, "A New Equation for the speed of sound in Natural
%            Waters", J. Acoust. Soc. Am. 56#4 (1974).
%
%          3) SOUNDSPEED(S,T,P,'chen') returns the sound speed (m/sec) given
%            given vectors of salinity (ppt), temperature (deg C), and 
%            PRESSURE (dbar) using the Chen and Millero equation:
%
%            Chen and Millero, "The Sound Speed in Seawater", J. Acoust. 
%            Soc. Am. 62 (1977), 1129-1135
%
%          4) SOUNDSPEED(S,T,P,'state') returns the sound speed (m/sec) given
%            given vectors of salinity (ppt), temperature (deg C), and 
%            PRESSURE (dbar) by using derivatives of the EOS80 equation of
%            state for seawater and the adiabatic lapse rate.
% 
%

%Notes: RP (WHOI) 3/dec/91
%        Added state equation ss

if (nargin<4), equation='mackenzie'; end

if (equation(1:3)=='mac')

% Fixed a small typo with the salinity stuff in the T*(S-35) term
% --RP 15/11/91

     c= 1.44896e3;     t= 4.591e0; 
    t2=-5.304e-2;     t3= 2.374e-4;
     s= 1.340e0;       d= 1.630e-2;
    d2= 1.675e-7;     ts=-1.025e-2;
   td3=-7.139e-13;  

   ssp=c+t*T+t2*T.*T+t3*T.*T.*T+s*(S-35.0)+d*D+d2*D.*D+ts*T.*(S-35.0) ...
        +td3*T.*D.*D.*D;

elseif (equation(1:3)=='del')

% Del grosso uses pressure in kg/cm^2. To get to this from dbars we must 
% divide by "g". From the UNESCO algorithms (referring to ANON (1970) 
% BULLETIN GEODESIQUE) we have this formula for g as a function of latitude 
% and pressure. We set latitude to 45 degrees for convenience!
      XX=sin(45*pi/180);
      GR = 9.780318*(1.0+(5.2788E-3+2.36E-5*XX)*XX) + 1.092E-6*D;
      P=D./GR;
% This is from VSOUND.f:
    
      C000 = 1402.392;
      DCT = (0.501109398873e1-(0.550946843172e-1 - 0.221535969240e-3*T).*T).*T;
      DCS = (0.132952290781e1 + 0.128955756844e-3*S).*S;
      DCP = (0.156059257041e0 + (0.244998688441e-4 - 0.883392332513e-8*P).*P).*P;
      DCSTP = -0.127562783426e-1*T.*S + 0.635191613389e-2*T.*P +0.265484716608e-7*T.*T.*P.*P ...
 - 0.159349479045e-5*T.*P.*P+0.522116437235e-9*T.*P.*P.*P - 0.438031096213e-6*T.*T.*T.*P;
      DCSTP=DCSTP - 0.161674495909e-8*S.*S.*P.*P + 0.968403156410e-4*T.*T.*S+ ...
   0.485639620015e-5*T.*S.*S.*P - 0.340597039004e-3*T.*S.*P;
   ssp= C000 + DCT + DCS + DCP + DCSTP;


elseif (equation(1:3)=='che')
      P0=D;
% This is copied directly from the UNESCO algorithmms, with some minor changes (like adding
% ";" and changing "*" to ".*") for Matlab.

% CHECKVALUE: SVEL=1731.995 M/S, S=40 (IPSS-78),T=40 DEG C,P=10000 DBAR

%   SCALE PRESSURE TO BARS
      P=P0/10.;
%**************************
      SR = sqrt(abs(S));
% S**2 TERM
      D = 1.727E-3 - 7.9836E-6*P;
% S**3/2 TERM
      B1 = 7.3637E-5 +1.7945E-7*T;
      B0 = -1.922E-2 -4.42E-5*T;
      B = B0 + B1.*P;
% S**1 TERM
      A3 = (-3.389E-13*T+6.649E-12).*T+1.100E-10;
      A2 = ((7.988E-12*T-1.6002E-10).*T+9.1041E-9).*T-3.9064E-7;
      A1 = (((-2.0122E-10*T+1.0507E-8).*T-6.4885E-8).*T-1.2580E-5).*T+9.4742E-5;
      A0 = (((-3.21E-8*T+2.006E-6).*T+7.164E-5).*T-1.262E-2).*T+1.389;
      A = ((A3.*P+A2).*P+A1).*P+A0;
% S**0 TERM
      C3 = (-2.3643E-12*T+3.8504E-10).*T-9.7729E-9;
      C2 = (((1.0405E-12*T-2.5335E-10).*T+2.5974E-8).*T-1.7107E-6).*T +3.1260E-5;
      C1 = (((-6.1185E-10*T+1.3621E-7).*T-8.1788E-6).*T+6.8982E-4).*T +0.153563;
      C0 = ((((3.1464E-9*T-1.47800E-6).*T+3.3420E-4).*T-5.80852E-2).*T+5.03711).*T+1402.388;
      C = ((C3.*P+C2).*P+C1).*P+C0;
% SOUND SPEED RETURN
      ssp = C + (A+B.*SR+D.*S).*S;

elseif(equation(1:3)=='sta'),
     P=D;
% (Copied somewhat from program EOSSPEED.F)
     [svan,sigma]=swstate(S,T,P);
     VOL=(1.)./(1000.+sigma);
%     DV/DP|ADIA = (DV/DP) AT CONSTANT T + ADIA.LAPSE RATE *
%                  (DV/DT) AT CONSTANT P
%     Note: factor of 10 is convert pressure from dB to Bars
     dVdP=swstate(S,T,P,'dP');
     dVdT=swstate(S,T,P,'dT');
     dVdPad=(dVdP+adiabattempgrad(S,T,P).*dVdT)*10;
%     C = V * SQRT ( 1/DV/DP| ADIA)
     ssp=VOL.*sqrt(abs( (1.e5)./dVdPad ));

else
   error('soundspeed: Unrecognizable equation specified!');
end                 
