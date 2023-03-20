function [ x ] = bisect_by_count( xL, xR )

% Returns all the root in an interval (xL, xR)
% Does all the roots in parallel
% [ xL, xR] does not need to isolate the roots
% This version works purely off the mode count

% note that every 10 bisections gives 3 decimal digits
% (10 bisections reduces the interval by 2^10 = 1024 which is 3 digits)

% Mike Porter, 11/2009

MaxBisections = 50;   % maximum number of bisections

[ ~, ModeCount1 ] = pekeris( xL );
[ ~, ModeCount2 ] = pekeris( xR );

Counter = ModeCount1 : -1 : ModeCount2 + 1;
M = length( Counter );

% Obtain initial guesses for x1 and x2
x1( 1 : M ) = xL;
x2( 1 : M ) = xR;

% Begin bisection
for J = 1 : MaxBisections
   xmid = x1 + ( x2 - x1 ) / 2;   % midpoint of bracketting interval
   [ ~, ModeCount ] = pekeris( xmid );   % evaluate the characteristic fn

   % replace left bracket
   ii = find( ModeCount >= Counter );
   x1( ii ) = xmid( ii );
   
   % replace right bracket
   ii = find( ModeCount < Counter );
   x2( ii ) = xmid( ii );
   
end

x = xmid;

return
