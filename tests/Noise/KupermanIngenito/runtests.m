% compares noise from the modal and spectral integral formulas
% Laurel Henderson and Mike Porter
% 11/2017

% The curves derived from the 'analytic' formula should agree closely with
% those computed from the field_noise
%
% In general KRAKEN (as opposed to KRAKENC) shows issues due to its perturbation theory
% In general the 'diag' option shows issues due to that approximation

% Differences in DownwdRef between KRAKENC and SCOOTER at 400 Hz are
% probably due to missing modes in KRAKENC

% The deviation between KRAKENC and SCOOTER at higher frequencies (sduct case) is
% probably due to the perturbation treatment of surface loss in KRAKENC
% The '.' option for Sduct caused major problems in KRAKENC

% Check sigma = 0.5 m in Sduct
% I've temporarily edited out the surface roughness in sduct ...
% That probably causes the NL taken directly from the shdfil to be too low
% (range disc of 100 km not big enough)
% Getting failure to converge in inverse iteration for a couple of modes in
% this case

% loop over all the test cases

filelist = [ ...
    'Iso      '
    'DownwdRef'
    'Sduct    ' ];

for icase = 1 : 3
    figure

    filename = deblank( filelist( icase, : ) )
    kraken(  filename  )

    filenameC = [ filename '_C' ]
    krakenc( filenameC )

    filenameS = [ deblank( filelist( icase, : ) ) '_S' ]
    GreenFile = [ filenameS '.grn' ];
    scooter( filenameS )

    filenameSnoise = [ filenameS 'noise' ]
    GreenFilenoise = [ filenameSnoise '.grn' ];
    scooter_nofield( filenameSnoise )

    for ifreq = 1 : 3
        freq = 2 ^ ifreq * 100;   % frequency in Hz

        % plot_noise

        Rmax_km = 1e9;   % max range of integration of noise sources

        ModeFile  = [ filename '.mod' ];
        ModeFileC = [ filenameC '.mod' ];

        sd = 0.5;  % depth of noise sources
        rd = 0 : 0.1 : 50;

        rho_SL_dB = 0;

        Component = 'P';

        %%
        % Spectral noise vs. depth
        % Note that SCOOTER has to be run with stabilizing attenuation disabled
        % ( TopOpt( 6 : 6 ) = '0' )

        NL = spectral_noise( GreenFilenoise, rho_SL_dB, sd, rd, freq );

        plot( NL, rd, 'k', 'LineWidth', ifreq );
        ll = legend( 'SCOOTER analytic' );

        % ll = legend( 'NM Full', 'NMC Full', 'NM Diag', 'NMC Diag' );
        set(ll,'Location','southeast')
        grid
        drawnow
        hold on

        % SCOOTER noise from shd

        C  = field_noise( [ filenameS '.shd.mat' ], sd, rd, freq );
        NL = 10 * log10( diag( C ) ) + rho_SL_dB;

        plot( NL, rd, '--k', 'LineWidth', ifreq );
        ll = legend( 'SCOOTER analytic', 'SCOOTER field' );
        set( ll, 'Location', 'southeast' )
        grid
        axis( [ 4 18 0 50 ] )
        drawnow

        %%
        % Modal noise vs. depth using full matrix
        NLC = modal_noise_full( ModeFileC, rho_SL_dB, sd, rd, freq, Rmax_km, Component );

        plot( NLC( :, end ), rd, 'b', 'LineWidth', ifreq );
        ll = legend( 'SCOOTER analytic', 'SCOOTER field', 'KRAKENC analytic' );
        ylabel( 'Receiver depth (m)' )
        xlabel( 'NL (dB)' )
        set(gca,'YDir','reverse')
        %tt = title(filename);
        tt = title( { deblank( filename ); [ 'Sd = ' num2str( sd ) ' m' ] } );
        set(tt,'Interpreter','none')

        % Modal noise vs. depth using diagonal terms only
        %       NL  = modal_noise_diag( ModeFile,  rho_SL_dB, sd, rd, Rmax_km, Component );
        %       NLC = modal_noise_diag( ModeFileC, rho_SL_dB, sd, rd, Rmax_km, Component );
        %       plot( NL(  :, end ), rd, 'g', 'LineWidth', 3 );
        %       plot( NLC( :, end ), rd, 'c', 'LineWidth', 3 );
        %
        % KRAKENC noise from shd

        C  = field_noise( [ filenameC '.shd.mat' ], sd, rd, freq );
        NL = 10 * log10( diag( C ) ) + rho_SL_dB;

        plot( NL, rd, '--b', 'LineWidth', ifreq );
        ll = legend( 'SCOOTER analytic', 'SCOOTER field', 'KRAKENC analytic', 'KRAKENC field' );
        set( ll, 'Location', 'southeast' )
        grid
        axis( [ 4 18 0 50 ] )
        drawnow
        %%
        if ( icase == 1 )
        % Modal noise vs. depth using full matrix
        NL  = modal_noise_full( ModeFile,  rho_SL_dB, sd, rd, freq, Rmax_km, Component );
        % we just take the value from the matrix NL for the largest range

        plot( NL(  :, end ), rd, 'c', 'LineWidth', ifreq );
        hold on

        % KRAKEN noise from shd

        C  = field_noise( [ filename '.shd.mat' ], sd, rd, freq );
        NL = 10 * log10( diag( C ) ) + rho_SL_dB;

        plot( NL, rd, '--c', 'LineWidth', ifreq );
        ll = legend( 'SCOOTER analytic', 'SCOOTER field', 'KRAKENC analytic', 'KRAKENC field', 'KRAKEN analytic', 'KRAKEN field' );
        set( ll, 'Location', 'southeast' )
        grid
        axis( [ 4 18 0 50 ] )
        drawnow
        end
    end   % next frequency

    % set( gca, 'Position', [ 2    2                       14.0       7.0 ] )
    set( gcf, 'Units', 'centimeters' )
    set( gcf, 'Position', [ 3 15 19.0 11.0 ] )
    print( filename, '-dpng' )
end   % next icase