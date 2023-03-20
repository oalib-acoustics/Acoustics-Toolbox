function map = haxby(m)
%
%HAXBY  Haxby color map
%   HAXBY(M) returns an M-by-3 matrix containing a colormap with Haxby's
%   colors, commonly used for displaying bathymetry data.
%   HAXBY, by itself, is the same length as the current colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(haxby)
%
%   Use
%             colormap(flipud(haxby))
%
%   for bathymetry data (positive downward).
%
%   Colormap is based on the colors used by W. F. Haxby's Gravity
%   field of World's oceans, 1985, developed for geoid and gravity maps.
%   The version used here is formed from a linear interpolation of
%   the GMT color table used by MB-System by David W. Caress and Dale N. Chayes.
%   <http://www.ldeo.columbia.edu/res/pi/MB-System>
%
%   See also HSV, GRAY, PINK, COOL, BONE, COPPER, FLAG, HOT
%   COLORMAP, RGBPLOT.

% Kelsey Jordahl
% Marymount Manhattan College
% Time-stamp: <Fri Oct 30 12:45:12 EDT 2009>

if nargin < 1, m = size(get(gcf,'colormap'),1); end

% mbm_grdplot Haxby color pallette
% ncolors=11;
% c=[ 37    57   175;    40   127   251;    50   190   255;   106   235   255;
%     138   236   174;   205   255   162;   240   236   121;   255   189    87;
%     255   161    68;   255   186   133;   255   255   255];

ncolors = 18;   % Gebco-blue-brown.cpt (just blue part)
c=[153     0       255   ;   
  136     17      255   ;   
  119     34      255   ;   
  102     51      255   ;   
  85      68      255   ;   
  68      85      255   ;   
  51      102     255   ;   
  34      119     255   ;   
  17      136     255   ;   
  0       153     255   ;   
  27      164     255   ;   
  54      175     255   ;   
  81      186     255   ;   
  108     197     255   ;   
  134     208     255   ;   
  161     219     255   ;   
  188     230     255   ;
  255     255     255  ] ;
%  200     235     255  ];   
%  210     241     255 ] ; 

% ncolors=17;                 %alternate haxby colormap.
% c= [10    0    121;
% 40    0    150;
% 20    5    175;
% 0    10    200;
% 0    25    212;
% 0    40    224;
% 26    102    240;
% 13    129    248;
% 25    175    255;
% 50    190    255;
% 68    202    255;
% 97    205    240;
% 106    215    225;
% 124    215    200;
% 138    226    174;
% 182    235    186;
% 255    255    255];

pp = 1 : (m-1)/(ncolors-1) : m;
r  = interp1( pp, c(:,1), 1:m );
g  = interp1( pp, c(:,2), 1:m );
b  = interp1( pp, c(:,3), 1:m );

map = [r' g' b']/255;

% make sure roundoff doesn't cause map values outside [ 0, 1 ]
map( map < 0 ) = 0;
map( map > 1 ) = 1;
