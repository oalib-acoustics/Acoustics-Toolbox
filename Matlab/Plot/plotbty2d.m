function plotbty2d( xyzfile )

% plot the bathymetry from an xyzfile

[ Bathy ] = LoadBathymetry( xyzfile );

% axes( handles.bty )
% set( handles.bty,'FontSize',14)
% 
% colorbar( 'delete' )
% axes( handles.TL )
% colorbar( 'delete' )
      
% Plot the bathymetry

% Turn 180 to 360 degrees East into -180 to -0 degrees (West)
if (Bathy.Lon > 180.)
  Bathy.Lon = Bathy.Lon - 360.;
end 

Bathy.depth( isnan( Bathy.depth )  ) = 0;

axis( [ min( Bathy.Lon ), max( Bathy.Lon ) min( Bathy.Lat ), max( Bathy.Lat ) ] )
% handles.btymap = pcolor( Bathy.Lon, Bathy.Lat,  Bathy.depth );
% handles.btymap = surfc(  Bathy.Lon, Bathy.Lat, -Bathy.depth );
handles.btymap = contourf( Bathy.Lon, Bathy.Lat, -Bathy.depth, -[ 0 : 250 : 6000 ] );

% shading interp
xlabel( 'Longitude' )
ylabel( 'Latitude' )
colormap( haxby( 28 ) )
caxis( [ -7000 0 ] );
colorbar

% % bty colorbar on the left
% cbarBTY = colorbar( 'Location', 'WestOutside' );
% 
% %%%%%%%%%%%%
% % Colorbar not crowding the y-lable
% % Modify Colorbar to a manual setting
% set(cbarBTY,'location','manual','ActivePositionProperty','OuterPosition')
% 
% % Reset the outerposition of the colorbar to be normalized from figure
% cbarBTYPos = get(cbarBTY,'OuterPosition');
% set(cbarBTY,'OuterPosition',[cbarBTYPos(3)/2 0 cbarBTYPos(3) 1])
% 
% %-------------------cpm 14/10/09
% % changed line below
% % set( cbarBTY, 'YLim', CLim1 );   %restrict colorbar axes limits
% 
% %  reverse y direction to have depth increasing downward
% set( cbarBTY, 'YLim', CLim1,'YDir','reverse' ); %restrict colorbar axes limits
% % and add a title to the colorbar
% set( get(cbarBTY,'Title'),'String','Depth (m)');
% set( get(cbarBTY,'Title'),'FontSize',14);
% %set(cbarBTY,'FontSize',[12]);                   % make TickLabels large, too
% %------------
% 
