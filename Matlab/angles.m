
function angles = angles( c1, c2 )

% ANGLES computes the take-off angle which is turned at sound speed c2.
%
% usage: angles( c1, c2 )
%
% Prompts for two sound speed values and uses Snell's law to calculate
% the launch angle (in degrees) which which would produce a horizontal ray at the depth of the second sound speed
% mbp 26 Feb. 2001

angles = acosd( c1 / c2 );

disp( 'Angle in degrees')
disp( angles )
