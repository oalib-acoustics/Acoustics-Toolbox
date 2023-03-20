function plottrc( trcfil )

% plot the Top Reflection Coefficient file used by Bellhop
% usage: plottrc( trcfil )
%
% MBP May 2002

global thetaTop RTop phiTop

% Read the reflection coefficient
brcfil = ' ';
readrc( trcfil, brcfil, 'F', ' ' )

[ AX, ~, H2 ] = plotyy( thetaTop, RTop, thetaTop, 180.0 / pi * phiTop, 'plot' );

set( get( AX( 1 ),'Ylabel'), 'String','|R|') 
set( get( AX( 2 ),'Ylabel'), 'String','angle (degrees)' )

xlabel( 'angle (degrees)' )
ylabel( 'Reflection Loss (dB)' )
title( 'Top Reflection Coefficient' )
