
% plot a series of snapshots of the acoustic field
% using a linear scale

for ifreq = 1:nfreq
    % open the file
    fid = fopen('feif_top.asc','r');

    % read

    pltitl = fgetl( fid );
    junk   = fgetl( fid );
    junk   = fscanf( fid, '%f', 1 );
    nsd    = fscanf( fid, '%i', 1 );
    nrd    = fscanf( fid, '%i', 1 );
    nrr    = fscanf( fid, '%i', 1 );

    sd     = fscanf( fid, '%f', nsd );
    rd     = fscanf( fid, '%f', nrd );
    rr     = fscanf( fid, '%f', nrr );

    isd = 11;
    for i = 1:isd
        temp1   = fscanf( fid, '%f', [ 2 * nrr, nrd ] );
        i, size(temp1)
    end

    isd = 10;
    for i = 1:isd
        temp2   = fscanf( fid, '%f', [ 2 * nrr, nrd ] );
    end
    for i = 1:isd
        temp3   = fscanf( fid, '%f', [ 2 * nrr, nrd ] );
    end
    for i = 1:isd
        temp4   = fscanf( fid, '%f', [ 2 * nrr, nrd ] );
    end

    zt = rd;
    taker = 1:nrr;
    rt = rr( taker );

    figure

    tlt = temp1( 1:2:2*nrr, : )';
    subplot( 4, 1, 1 )
    pcolor( rt, zt, tlt ); ...
        shading flat; colormap( gray ); colorbar; view( 0, -90 );
    hold on;
    hidden off;
    contour( rt, zt, tlt, 10 );
    axis('off');
    ylabel( 'Depth (m)' );
    title('t = 0.02 s')

    tlt = temp2( 1:2:2*nrr, : )';
    subplot( 4, 1, 2 )
    pcolor( rt, zt, tlt ); ...
        shading flat; colormap( gray ); colorbar; view( 0, -90 );
    hold on;
    hidden off;
    contour( rt, zt, tlt, 10 );
    axis('off');
    ylabel( 'Depth (m)' );
    title('t = 0.04 s')

    tlt = temp3( 1:2:2*nrr, : )';
    subplot( 4, 1, 3 )
    pcolor( rt, zt, tlt ); ...
        shading flat; colormap( gray ); colorbar; view( 0, -90 );
    hold on;
    hidden off;
    contour( rt, zt, tlt, 10 );
    axis('off');
    ylabel( 'Depth (m)' );
    title('t = 0.06 s')

    tlt = temp4( 1:2:2*nrr, : )';
    subplot( 4, 1, 4 )
    pcolor( rt, zt, tlt ); ...
        shading flat; colormap( gray ); colorbar; view( 0, -90 );
    hold on;
    hidden off;
    contour( rt, zt, tlt, 10 );
    xlabel( 'Range (m)' ); ylabel( 'Depth (m)' );
    title('t = 0.08 s')


    set(1,'PaperPosition', [ 0.25 0.00 6.5 7.0 ] )
    print -deps sparc.ps
end