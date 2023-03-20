function [ x ] = bisect_final( xL, xR )

% Returns the root using an isolating interval (xL, xR) for each root
% Does all the roots in parallel
% This version uses the sign of the function, as well as the mode count

% Mike Porter, 11/2009

MaxBisections = 50;   % maximum number of bisections

[ ~, ModeCount ] = pekeris( xR( 1 ) );

% Obtain initial guesses for x1 and x2
x2 = xR;
x1 = xL;

[ f1, ModeCount ] = pekeris( x1 );
[ f2, ModeCount ] = pekeris( x2 );

% Begin bisection (allowing no more than MaxBIS bisections per mode)
for J = 1 : MaxBisections
   xmid = x1 + ( x2 - x1 ) / 2;
   [ fmid, ModeCount ] = pekeris( xmid );
   
   %[ sign( f1 ) sign( f2 ) sign( fmid ) ];
   % replace left bracket
   ii = find( sign( fmid ) == sign( f1 ) );
   x1( ii ) = xmid( ii );
   f1( ii ) = fmid( ii );
   
   % replace right bracket
   ii = find( sign( fmid ) == sign( f2 ) );
   x2( ii ) = xmid( ii );
   f2( ii ) = fmid( ii );
   
end

x = xmid;

return
