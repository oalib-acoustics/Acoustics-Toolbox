function [ z, psiScaled ] = normiz( psi, X, B3, B4, rho, Mater )

% Convert the eigenvector to (U, W, TAUZX, TAUZZ) and normalize
% density is assumed constant in acoustic layers
%
% Halfspaces are not allowed in the Matlab version of Krakel since it
% requires a conventional algebraic eigenvalue problem for the Matlab
% solver
%
% The vector z must be contained in the acoustic layers because the
% interpolation routine does not deal with the four independent components of
% the stress-displacement vector for an elastic medium.

% As written, the code returns everything in psiScaled
% Separately, it extracts u, v, tauzx, tauzz, p but does not return those values

global SSP
global omega
global Loc h

ik     = 1i * sqrt( X );
NMat   = length( psi );

omega2 = omega^2;

% Scale down the mode

SqNorm2 = norm( psi, 1 );
psi     = psi / SqNorm2;

% Loop to compute norm of eigenvector

SqNorm = 0.0;
jj     = 1;
kk     = 1;
psiScaled   = zeros( size( psi ) );

z( 1 ) = 0.0;    % assumes DepthT = 0 for now

for Medium = 1 : SSP.NMedia
    rhoM = rho( Loc( Medium ) + 1 );
    for ii = 1 : SSP.N( Medium ) + 1
        if ( strcmp( Mater( Medium, : ), 'ELASTIC ' ) ) % elastic layer
            R1 = psi( kk     );
            R2 = psi( kk + 1 );
            R3 = psi( kk + 2 );
            R4 = psi( kk + 3 );

            Contrib = h( Medium ) * ( -X * B3( jj ) * R1^2 + B4( jj ) * R1 * R4 - R2 * R3 );

            % Load components sequentially into psiS
            psiScaled( kk     ) = ik * psi( kk     );
            psiScaled( kk + 1 ) =      psi( kk + 1 );
            psiScaled( kk + 2 ) = ik * psi( kk + 2 );
            psiScaled( kk + 3 ) =      psi( kk + 3 );

            % Load components into separate vectors
            u(     ii, Medium ) = ik * psi( kk     );
            v(     ii, Medium ) =      psi( kk + 1 );
            tauzx( ii, Medium ) = ik * psi( kk + 2 );
            tauzz( ii, Medium ) =      psi( kk + 3 );

            kk = kk + 4;
        else  % acoustic layer
            Contrib         = -h( Medium ) * psi( kk )^2 / ( rhoM * omega2 );
            psiScaled( kk ) = -psi( kk );
            p( ii, Medium ) = -psi( kk );

            kk = kk + 1;
            z( kk ) = z( kk - 1 ) + h( Medium );   % position for next point
        end

        if ( ii == 1 || ii == SSP.N( Medium ) + 1 )
            Contrib = 0.5 * Contrib;
        end
        SqNorm = SqNorm + Contrib;
        jj = jj + 1;
    end
end

z = z( 1 : end - 1 );

% Scale the mode

RN       = -omega2 * SqNorm;
ScaleFac = 1.0 / sqrt( RN );
psiScaled( 1 : NMat ) = ScaleFac * psiScaled( 1 : NMat );
