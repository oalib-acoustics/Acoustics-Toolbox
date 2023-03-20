function [ w, Ix ] = weight( x, xTab )

% Computes weights for linear interpolation
%     Given 
%        x(*)    abscissas
%        xTab(*) points for tabulation
%
%     Compute
%        w(*)    weights for linear interpolation
%        Ix(*)   indices for    "         "
%

validateattributes( x, { 'numeric' }, { 'vector', 'increasing' } )


% pre-allocate output vectors for efficiency
Nx    = length( x    );
NxTab = length( xTab );
Ix    = zeros( NxTab, 1 );
w     = zeros( NxTab, 1 );

% Quick return if just one X value for interpolation ***

if ( Nx == 1 )
    w(  1 ) = 0.0;
    Ix( 1 ) = 1;
    return
end

L = 1;

for IxTab = 1 : NxTab   % Loop over each point for which the weights are needed
    
    % search for index, L, such that [X(L), X(L+1)] brackets rcvr depth
    % (this could also be done as min( xTab - x )
    while ( xTab( IxTab ) > x( L + 1 ) && L < Nx - 1 )
        L = L + 1;
    end
    Ix( IxTab ) = L;                                                    % here's the output index
    w(  IxTab ) = ( xTab( IxTab ) - x( L ) ) / ( x( L+1 ) - x( L ) );   % here's the output weight
end