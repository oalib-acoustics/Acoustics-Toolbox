function march( SSP, c2R, c2I, rho, x, Ik, deltak, rr )

% Solve system for a sequence of k-values

global H
global Pulse fMax tStart alpha beta V
global tout deltat
global cLow cHigh
global omega Bdry
global WS isz WR irz Z
global Green RTSrd RTSrr

MaxIT   = 1e6;
deltat2 = deltat^2;

% March a single spectral component forward in time
NTot1 = sum( SSP.N ) + 1;
AD0   = zeros( NTot1, 1 );
AE0   = zeros( NTot1, 1 );
AD1   = zeros( NTot1, 1 );
AE1   = zeros( NTot1, 1 );
AD2   = zeros( NTot1, 1 );
AE2   = zeros( NTot1, 1 );
Ke    = zeros( 2, 2 );
Ce    = zeros( 2, 2 );
Me    = zeros( 2, 2 );

% Assemble mass and stiffness matrices

rkT  = sqrt( x );

Node = 1;
L    = 1;

for medium = 1 : SSP.NMedia
    hElt = H( medium );

    for I = 1 : SSP.N( medium )
        rhoElt = ( rho( L ) + rho( L + 1 ) ) / 2.0;
        c2RElt = ( c2R( L ) + c2R( L + 1 ) ) / 2.0;
        c2IElt = ( c2I( L ) + c2I( L + 1 ) ) / 2.0;
        rhoh   = rhoElt * hElt;

        % Elemental stiffness matrix ( +k2 term )
        y = 1.0 + 1i * rkT * V * c2IElt;
        z = hElt * x * ( y - V * V / c2RElt ) / rhoElt / 6.0;
        Ke( 1, 1 ) =  y / rhoh + (3.0-alpha) * z;
        Ke( 1, 2 ) = -y / rhoh +      alpha  * z;
        Ke( 2, 2 ) =  y / rhoh + (3.0-alpha) * z;

        % Elemental damping matrix
        xx = hElt * x / rhoElt / 6.0;
        z  = hElt * 2.0 * 1i * rkT * V / c2RElt / rhoElt / 6.0;
        Ce( 1, 1 ) = c2IElt * ( 1.0/rhoh + (3.0-alpha) * xx) + (3.0-alpha) * z;
        Ce( 1, 2 ) = c2IElt * (-1.0/rhoh +      alpha  * xx) +      alpha  * z;
        Ce( 2, 2 ) = c2IElt * ( 1.0/rhoh + (3.0-alpha) * xx) + (3.0-alpha) * z;

        % Elemental mass matrix
        Me( 1, 1 ) = (3.0-alpha) * hElt / ( rhoElt * c2RElt ) / 6.0;
        Me( 1, 2 ) =      alpha  * hElt / ( rhoElt * c2RElt ) / 6.0;
        Me( 2, 2 ) = (3.0-alpha) * hElt / ( rhoElt * c2RElt ) / 6.0;

        % A2 matrix
        AD2( Node   ) = AD2( Node ) + ...
            Me(1,1) + 0.5 * deltat * Ce(1,1) + beta * deltat2 * Ke(1,1);
        AE2( Node+1 ) = Me(1,2) + 0.5 * deltat * Ce(1,2) + beta * deltat2 * Ke(1,2);
        AD2( Node+1 ) = Me(2,2) + 0.5 * deltat * Ce(2,2) + beta * deltat2 * Ke(2,2);

        % A1 matrix
        AD1( Node   ) = AD1( Node ) + ...
            2.0*Me(1,1) - ( 1.0 - 2.0*beta ) * deltat2 * Ke(1,1);
        AE1( Node+1 ) = 2.0*Me(1,2) - ( 1.0 - 2.0*beta ) * deltat2 * Ke(1,2);
        AD1( Node+1 ) = 2.0*Me(2,2) - ( 1.0 - 2.0*beta ) * deltat2 * Ke(2,2);

        % A0 matrix
        AD0( Node   ) = AD0( Node ) ...
            -Me(1,1) + 0.5 * deltat * Ce(1,1) - beta * deltat2 * Ke(1,1);
        AE0( Node+1 ) = -Me(1,2) + 0.5 * deltat * Ce(1,2) - beta * deltat2 * Ke(1,2);
        AD0( Node+1 ) = -Me(2,2) + 0.5 * deltat * Ce(2,2) - beta * deltat2 * Ke(2,2);

        Node = Node + 1;
        L    = L    + 1;
    end
    L = L + 1;
end

%factor( NTot1, AD2, AE2, RV1, RV2, RV3, RV4 )  % * Factor A2 *

% * Initialize pressure vectors *
U0 = zeros( NTot1, 1 );
U1 = zeros( NTot1, 1 );

% Initializate parameters for bandpass filtering of the source time series
if ( Pulse(4:4) == 'L' || Pulse(4:4) == 'B' )
    fLoCut = rkT * cLow / ( 2.0 * pi );
else
    fLoCut = 0.0;
end

if ( Pulse(4:4) == 'H' || Pulse(4:4) == 'B' )
    fHiCut = rkT * cHigh / ( 2.0 * pi );
else
    fHiCut = 10.0 * fMax;
end

time    = 0.0;
IniFlag = 1;   % need to tell source routine to refilter for each new value of k

% Forward, march
n  = length( AE0 );
AF0 = [ AE0( 2 : end ); 0 ];    % to conform to Matlab convention relating diagonals to the matrix
AF1 = [ AE1( 2 : end ); 0 ];
A0 = spdiags( [ AF0 AD0 AE0 ], -1:1, n, n );
A1 = spdiags( [ AF1 AD1 AE1 ], -1:1, n, n );
itout = 1;

for Itime = 1 : MaxIT
    time  = tStart + ( Itime - 1 ) * deltat;
    NTot1 = length( U0 );

    % ** Form U2TEMP = A1*U1 + A0*U0 + S **
    U2 = A1 * U1 + A0 * U0;

    %Ntpts = 1;
    %timeV( 1 ) = time; 	% source expects a vector, not a scalar
    %source( timeV, ST, sd, Nsd, Ntpts, omega, fLoCut, fHiCut, Pulse, PulseTitle, IniFlag );
    [ ST, PulseTitle ] = cans( time, omega, Pulse );
    
    ST = deltat2 * ST * exp( -1i * rkT * V * time );
    js = isz;    % assumes source in first medium
    U2( js   ) = U2( js   ) + ( 1.0 - WS ) * ST;
    U2( js+1 ) = U2( js+1 ) +         WS   * ST;

    % Solve A2*U2 = U2TEMP
    if ( alpha == 0.0 && beta == 0.0 )   % explicit solver
        U2 = U2 ./ AD2;
    else
        backsb( NTot1, RV1, RV2, RV3, RV4, U2 );   % implicit solver
    end
    
    % Do a roll
    U0 = U1;
    U1 = U2;

    % Boundary conditions (natural is natural)
    if ( Bdry.Top.Opt(2:2) == 'V' )
        U1( 1 )     = 0.0;
    end
    if ( Bdry.Bot.Opt(1:1) == 'V' )
        U1( NTot1 ) = 0.0;
    end

    % Extract solution (snapshot, vertical or horizontal array)
    while ( ( itout <= length( tout ) ) && ( time + deltat >= tout( itout ) ) )
        % above can access tout( ntout + 1 ) which is outside array dimension but result is irrelevant
        wt = ( tout( itout ) - time ) / deltat;  % Weight for temporal interpolation:

        switch ( Bdry.Top.Opt(5:5) )
            case ( 'S' ) % snapshot
                const = exp( 1i * rkT * V * tout( itout ) );
                UT1 = U0( irz ) + WR .* ( U0( irz+1 ) - U0( irz ) );   % Linear interpolation in depth
                UT2 = U1( irz ) + WR .* ( U1( irz+1 ) - U1( irz ) );

                % note Green is 4-d: Green( Ntout, Nsd, Nrd, Nk )
                Green( itout, 1, :, Ik ) = const * ( UT1 + wt * ( UT2 - UT1 ) );   % Linear interpolation in time

            case ( 'D' ) % RTS (vertical array)
                const =  sqrt( 2.0 ) * deltak * sqrt( rkT ) * exp( 1i*( rkT * ( V * tout( itout ) - rr( 1 ) ) + pi / 4.0 ) );
                UT1 = U0( irz ) + WR .* ( U0( irz+1 ) - U0( irz ) );   % Linear interpolation in depth
                UT2 = U1( irz ) + WR .* ( U1( irz+1 ) - U1( irz ) );
                U   = UT1 + wt * ( UT2 - UT1 );                        % Linear interpolation in time
                RTSrd( :, itout ) = RTSrd( :, itout ) + real( const * U );

            case ( 'R' ) % RTS (horizontal array)
                if ( time >= tout( 1 ) )
                    % Linear interpolation in depth
                    I     = Ird( 1 );
                    T     = ( rd( 1 ) - Z( I ) ) / ( Z( I + 1 ) - Z( I ) );
                    UT1   = U0( I ) + T * ( U0( I+1 ) - U0( I ) );
                    UT2   = U1( I ) + T * ( U1( I+1 ) - U1( I ) );
                    U     = UT1 + wt * ( UT2 - UT1 );                       % Linear interpolation in time
                    const = sqrt( 2.0 ) * exp( 1i * rkT * V * time );
                    RTSrr( :, itout ) = RTSrr( :, itout ) + real( deltak * U *  ...
                        const * exp( 1i * ( -rkT * r + pi / 4.0 ) ) * sqrt( rkT / rr ) );
                end
        end
        itout = itout + 1;
    end
    if ( itout > length( tout ) )    % jump out of while loop and do next wavenumber
        break;
    end
end