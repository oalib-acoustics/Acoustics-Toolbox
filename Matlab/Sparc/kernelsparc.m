function kernelsparc( filename, SSP, c2R, c2I, rho, crosst, cMin, cMax, k, freq, Pos, rr, Bdry, PlotTitle )

% Solve system for a sequence of k-values

global tMult
global Green RTSrd RTSrr tout deltat

% pre-allocate for efficiency
Nrz    = length( Pos.r.z );
Nr     = length( rr );
Nk     = length( k );
Ntout  = length( tout );
deltak = ( k( end ) - k( 1 ) ) / Nk;

switch( Bdry.Top.Opt( 5 : 5 ) )
    case ( 'D')
        RTSrd = zeros( Nrz, Ntout );
    case ( 'R')
        RTSrr = zeros( Nr,  Ntout );
    otherwise
        Green = zeros( Ntout, 1, Nrz, Nk );   % Nsz=1 (matrix is compatible with SCOOTER)
end

deltat = tMult / sqrt( 1.0 / crosst^2 + ( 0.5 * cMax * k( Nk ) )^2 ); % Courant condition to set time step

fprintf( '\n\nTime step (based on CFL condition) = %f', deltat )
fprintf( '\n\nEstimated fl. pt. ops (millions) = %f \n', ( tout( end ) / deltat ) * Nk * length( c2R ) / 25000 )

% *** Loop over spectral components ***

for ik = 1 : Nk
    x = k( ik )^2;
    deltat = tMult / sqrt( 1.0 / crosst^2 + ( 0.5 * cMax * k( ik ) )^2 );
    if ( mod( ik, 10 ) == 0 )   % print every 10th wavenumber
       fprintf( '\n ik, Nk %i %i', ik, Nk )
    end
    march( SSP, c2R, c2I, rho, x, ik, deltak, rr );  % March that component for all time
end

freqVec( 1 ) = freq;

switch ( Bdry.Top.Opt( 5 : 5 ) )
    case ( 'S' )          % *** Case of a snapshot ***
        PlotType = 'rectilin  ';
        Pos.r.r = 2 * pi * freq ./ k.';   % store phasespeeds in slot usually used for receiver ranges
        freqVec = tout;  % store times       in slot usual used for frequencies

        % note pressure is 4-d: pressure( Ntout, Nsd, Nrd, Nk )
        pressure( :, :, :, : ) = Green;  % store Green's function in slot usually used for pressure
        freq0 = freq;
        atten = 0;
        save( [ filename '.grn.mat' ], 'PlotTitle', 'PlotType', 'freqVec', 'freq0', 'atten', 'Pos', 'pressure' )

    case ( 'D' )   % *** Write out RTS (vertical array) *
        Scale = 1.0 / sqrt( pi * rr( 1 ) );
        Pos.r.r = k;   % store wavenumbers is slot usually used for receiver ranges
        RTS = Scale * RTSrd;
        save( [ filename '.rts.mat' ], 'PlotTitle', 'freqVec', 'Pos', 'tout', 'RTS' )

    case ( 'R' )   %  *** Write out RTS (horizontal array) ***
        Pos.r.r = k;   % store wavenumbers is slot usually used for receiver ranges
        RTS = RTSrd;
        save( [ filename '.rts.mat' ], 'PlotTitle', 'freqVec', 'Pos', 'tout', 'RTS' )

end
