
% plot the rays associated with horizontally-refracted modes.
% MBP/Scripps/11/97

M = 1;	% number of modes propagating at source

fid = fopen( 'RAYFIL' );

%figure
xlabel( 'EW coordinate (km)')
ylabel( 'NS coordinate (km)')
zlabel( 'Depth (m)' )
%title( 'Key West Experiment: Ray paths, 50 Hz, Source at (-56, 30)' )
%grid on
hold on

for iray = 1 : 1000
  for mode = 1 : M
    Npts = fscanf( fid, '%i', 1 );
    temp = fscanf( fid, '%f', [ 2, Npts ] );

    xray = temp( 1, : ) / 1000;
    yray = temp( 2, : ) / 1000;
    col = 'krgbymc';	% colors to cycle through
    plot( xray, yray, col( mode ) )
  end		% next mode
  drawnow
end		% next ray

view( 2 )
view( 3 )
