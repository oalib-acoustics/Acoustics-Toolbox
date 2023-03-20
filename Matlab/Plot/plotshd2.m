
% plot a series of snapshots of the acoustic field
% using a linear scale

[ PlotTitle, PlotType, freq, atten, Pos, pressure ] = read_shd( filename );
tl = -squeeze( 20 * log10( abs( pressure ) ) );

zt = Pos.r.z;
rt = Pos.r.r;

figure
tej = flipud( colormap( jet ) );

tlt = squeeze( tl( 1, :, : ) );
subplot( 4, 1, 1 )
pcolor( rt, zt, tlt ); ...
shading flat; colormap( tej ); caxis( [ 40 80 ] ); colorbar; view( 0, -90 );
xlabel( 'Range (km)' );
ylabel( 'Depth (m)' );
title('Sd = 10 m')

tlt = squeeze( tl( 2, :, : ) );
subplot( 4, 1, 2 )
pcolor( rt, zt, tlt ); ...
shading flat; colormap( tej ); caxis( [ 40 80 ] ); colorbar; view( 0, -90 );
ylabel( 'Depth (m)' );
title('Sd = 20 m')

tlt = squeeze( tl( 3, :, : ) );
subplot( 4, 1, 3 )
pcolor( rt, zt, tlt ); ...
shading flat; colormap( tej ); caxis( [ 40 80 ] ); colorbar; view( 0, -90 );
ylabel( 'Depth (m)' );
title('Sd = 40 m')

tlt = squeeze( tl( 4, :, : ) );
subplot( 4, 1, 4 )
pcolor( rt, zt, tlt ); ...
shading flat; colormap( tej ); caxis( [ 40 80 ] ); colorbar; view( 0, -90 );
xlabel( 'Range (m)' ); ylabel( 'Depth (m)' );
title('Sd = 60 m')

%title( deblank( PlotTitle ) )

%set(1,'PaperPosition', [ 0.25 0.00 6.5 7.0 ] )
orient tall
print -dpng plotshd2.png
