function plotbrc( brcfil )

% plot the Bottom Reflection Coefficient file used by Bellhop
% usage: plotbrc( brcfil )
%
% MBP May 2002

global thetaBot RBot phiBot NBotPts thetaTop RTop phiTop NTopPts PlotTitle freqRC

% Read the reflection coefficient
trcfil = ' ';
readrc( trcfil, brcfil, ' ', 'F' )

%[ AX, ~, H2 ] = plotyy( thetaBot, 20 * log10( RBot ), thetaBot, 180.0 / pi * phiBot, 'plot' );
%[ AX, ~, H2 ] = plotyy( thetaBot,        abs( RBot ), thetaBot, 180.0 / pi * phiBot, 'plot' );
figure
plot( thetaBot, -20 * log10( RBot ), 'b', 'LineWidth', 3  );
xlabel( 'angle (degrees)' )
ylabel( 'Reflection Loss (dB)' )

figure
plot( thetaBot, abs( RBot ), 'b', 'LineWidth', 3  );
xlabel( 'angle (degrees)' )
ylabel( 'Reflection Loss' )

%set( get( AX( 1 ),'Ylabel'), 'String', '|R|') 
%set( get( AX( 2 ),'Ylabel'), 'String', 'angle (degrees)' )


%title( PlotTitle )
%title( { deblank( PlotTitle ); [ 'Freq = ' num2str( freqRC ) ' Hz' ] } )
