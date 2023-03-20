function [ Asparse, Bsparse ] = setAB(  B1, B2, B3, B4, rho, Mater )

% Initializes arrays defining the finite-difference equations

global Bdry SSP omega
global Loc h

X      = 0;
omega2 = omega^2;

% Compute size of matrix
NMat = 0;
for Medium = 1 : SSP.NMedia
    switch Mater( Medium, : )
        case 'ACOUSTIC'
            NMat = NMat + SSP.N( Medium ) + 1;
        case 'ELASTIC '
            NMat = NMat + 4 * ( SSP.N( Medium ) + 1 );
    end
end

A = zeros( 13, NMat ); % Zero out the matrix
B = zeros( 13, NMat ); % Zero out the matrix

% MAIN LOOP ( over each Medium )

J = 0;
for Medium = 1 : SSP.NMedia
    if ( Medium ~= 1 )

        % Acousto-elastic interface
        if ( strcmp( Mater( Medium - 1, : ), 'ACOUSTIC' ) && strcmp( Mater( Medium,: ), 'ELASTIC ' ) )
            % R2 Condition
            L = L + 1;
            J = J + 1;
            h_rho = h( Medium - 1 ) * rho( Loc( Medium - 1 ) + 1 );
            A( 10, J - 1 ) = 1.0;
            A(  9, J     ) = 0.5 * B1( L );
            B(  9, J     ) = 0.5 * h( Medium - 1 )^2;
            A(  7, J + 2 ) = h_rho * omega2;
            % R3 Condition
            A( 10, J     ) = 0.0;
            A(  7, J + 3 ) = 1.0;
            % R4 Condition
            A( 11, J     ) = 1.0;
            A(  7, J + 4 ) = 1.0;
        end

        % Elasto-acoustic interface

        if ( strcmp( Mater( Medium - 1, : ), 'ELASTIC ' ) && strcmp( Mater( Medium, : ), 'ACOUSTIC' ) )
            L = Loc( Medium ) + 1;
            J = J + 5;
            % R3 Condition
            A(9, J - 2) = 1.0;
            % R4 Condition
            A(9, J - 1) = 1.0;
            A(8, J    ) = 1.0;
            % R2 Condition
            h_rho = h( Medium ) * rho( Loc( Medium ) + 1 );
            A(12, J - 3 ) = -h_rho * omega2;
            A( 9, J     ) = 0.5 * B1( L );
            B( 9, J     ) = 0.5 * h( Medium )^2;
            A( 8, J + 1 ) = 1.0;
        end

        % Acoustic/acoustic interface

        if ( strcmp( Mater( Medium - 1, : ), 'ACOUSTIC' ) && strcmp( Mater( Medium, : ), 'ACOUSTIC' ) )
            % Continuity of normal displacement
            L = L + 1;
            J = J + 1;
            h_rho = h( Medium - 1 ) * rho( Loc( Medium - 1 ) + 1 );
            A( 10, J - 1 ) = 1.0 / h_rho;
            A(  9, J   ) = 0.5 * B1( L ) / h_rho;
            B(  9, J   ) = 0.5 * h( Medium - 1)^2 / h_rho;
            L = L + 1;
            h_rho = h( Medium ) * rho( Loc( Medium ) + 1 );
            A(  8, J + 1 ) = 0.5 * B1( L ) / h_rho;
            B(  8, J + 1 ) = 0.5 * h( Medium )^2 / h_rho;
            A(  7, J + 2 ) = 1.0 / h_rho;
            % Continuity of pressure
            J = J + 1;
            A( 10, J - 1 ) = 1.0;
            A(  9, J     ) = -1.0;
        end

        % Elastic/elastic interface

        if ( strcmp( Mater( Medium - 1, : ), 'ELASTIC ' ) && strcmp( Mater( Medium, : ), 'ELASTIC ' ) )
            % Continuity of everything
            J = J + 1;
            A( 11, J     ) =  1.0;
            A(  7, J + 4 ) = -1.0;
            J = J + 1;
            A( 11, J     ) =  1.0;
            A(  7, J + 4 ) = -1.0;
            J = J + 1;
            A( 11, J     ) =  1.0;
            A(  7, J + 4 ) = -1.0;
            J = J + 1;
            A( 11, J     ) =  1.0;
            A(  7, J + 4 ) = -1.0;
        end
    end

    % this section for the volume
    
    switch Mater( Medium,: )
        case 'ACOUSTIC' % Acoustic section
            L = Loc( Medium ) + 1;
            if ( Medium == 1 )
                J = 1;
            end

            for II = 1 : SSP.N( Medium ) - 1
                J = J + 1;
                L = L + 1;
                A( 10, J - 1 ) = 1.0;
                A(  9, J     ) = B1( L );
                B(  9, J     ) = h( Medium ) * h( Medium );
                A(  8, J + 1 ) = 1.0;
            end
        case 'ELASTIC '  % Elastic section

            L = Loc( Medium );
            two_by_h = 2.0 / h( Medium );
            for II = 1 : SSP.N( Medium )
                L = L + 1;
                J = J + 1;
                A( 11, J ) = two_by_h;
                B( 12, J ) = -B4( L );
                A( 13, J ) = -rho( L );
                B( 13, J ) = -B3( L );
                J = J + 1;
                A( 10, J ) = -1.0;
                A( 11, J ) = two_by_h;
                A( 13, J ) = -rho( L );
                J = J + 1;
                A(  9, J ) =  B1( L );
                A( 11, J ) = two_by_h;
                B( 12, J ) = -1;
                J = J + 1;
                A(  9, J ) =  B2( L );
                A( 10, J ) = -B4( L );
                A( 11, J ) = two_by_h;

                L = L + 1;
                J = J + 1;
                A( 7, J ) = -two_by_h;
                B( 8, J ) = -B4( L );
                A( 9, J ) = -rho( L );
                B( 9, J ) = -B3( L );
                J = J + 1;
                A( 6, J ) = -1.0;
                A( 7, J ) = -two_by_h;
                A( 9, J ) = -rho( L );
                J = J + 1;
                A( 5, J ) =  B1( L );
                A( 7, J ) = -two_by_h;
                B( 8, J ) = -1;
                J = J + 1;
                A( 5, J ) =  B2( L );
                A( 6, J ) = -B4( L );
                A( 7, J ) = -two_by_h;
                J = J - 4;
                L = L - 1;
            end
    end
end   % next Medium

% Top BC

switch Mater( 1,: )
    case 'ELASTIC '  % Elastic medium
        switch Bdry.Top.Opt( 1 : 1 )
            case 'R'  % Rigid top
                A( 9, 1 ) = 1.0;
                A( 9, 2 ) = 1.0;
            case 'V'  % Free (vacuum) top
                A( 7, 3 ) = 1.0;
                A( 7, 4 ) = 1.0;
        end
    case 'ACOUSTIC' % Acoustic medium
        % Top impedance
        [ F, G, IPow ] = bcimp( B1, B2, B3, B4, rho, X, 'TOP', Bdry.Top );
        if ( G == 0.0 )  % Free (vacuum) top
            A( 9, 1 ) = 1.0;
            A( 8, 2 ) = 0.0;
        else
            h_rho = rho(1) * h( 1 );
            A( 9, 1 ) = 0.5 * B1( 1 ) / h_rho + F / G;
            B( 9, 1 ) = 0.5 * h( 1 )^2 / h_rho - F / G;
            A( 8, 2 ) = 1.0 / h_rho;
        end
end

% Bottom BC

switch Mater( SSP.NMedia, : )
    case 'ELASTIC '  % Elastic medium
        switch Bdry.Bot.Opt(1:1 )
            case 'R'  % Rigid bottom
                A( 11, NMat - 3 ) = 1.0;
                A( 11, NMat - 2 ) = 1.0;
            case 'V'  % Free (vacuum) bottom
                A(  9, NMat - 1 ) = 1.0;
                A(  9, NMat     ) = 1.0;

            case 'A'  % Elastic bottom
                error( 'need to fix eigenvalue in this section' )
%                 GAMP2 = X - omega2 / CPB^2;
%                 GAMS2 = X - omega2 / CSB^2;
%                 GAMP = SQRT( GAMP2 );
%                 GAMS = SQRT( GAMS2 );
%                 RMU = rhoB * CSB^2;
%                 FAC = -RMU / ( GAMS * GAMP - X );
% 
%                 A( 11, NMat - 3 ) = FAC * GAMP * (X - GAMS2);
%                 A( 10, NMat - 2 ) = FAC * (2.0 * GAMS * GAMP - (GAMS2 + X));
%                 A( 12, NMat - 3 ) = X * A(10, NMat - 2);
%                 A( 11, NMat - 2 ) = FAC * GAMS * (X - GAMS2);
%                 A(  9, NMat - 1 ) = 1.0;
%                 A(  9, NMat   ) = 1.0;
        end
    case 'ACOUSTIC' % Acoustic medium
        % Bottom impedance
        [ F, G, IPow ] = bcimp( B1, B2, B3, B4, rho, X, 'BOT', Bdry.Bot );
        
        if ( G == 0.0 )  % Free (vacuum) bottom
            A( 10, NMat - 1 ) = 0.0;
            A(  9, NMat   ) = 1.0;
        else
            L = L + 1;
            h_rho = h( SSP.NMedia ) * rho( Loc( SSP.NMedia ) + 1 );
            A( 10, NMat - 1) = 1.0 / h_rho;
            A(  9, NMat    ) = 0.5  *  B1( L ) / h_rho - F / G;
            B(  9, NMat    ) = 0.5  *  h( SSP.NMedia )^2 / h_rho - F / G;
        end
end

Asparse = spdiags( A( 5:13, : ).', 4:-1:-4, NMat, NMat );
Bsparse = spdiags( B( 5:13, : ).', 4:-1:-4, NMat, NMat );
