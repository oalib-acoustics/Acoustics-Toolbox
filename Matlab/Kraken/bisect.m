function [ xL, xR ] = bisect( xMin, xMax )

% Returns an isolating interval (xL, xR) for each eigenvalue

% Mike Porter, 11/2009

MaxBisections = 50;

[ ~, ModeCount ] = pekeris( xMax );
NZero1 = ModeCount;   % mode number at xMax (usually 0)

[ ~, ModeCount ] = pekeris( xMin );
M = ModeCount - NZero1;   % number of modes

xL( 1 : M ) = xMin;   % initial left  boundary
xR( 1 : M ) = xMax;   % initial right boundary

if ( M == 1 )
   return   % quick exit if only one mode is sought
end

for Mode = 1 : M - 1  % loop over eigenvalue
   
   % Obtain initial guesses for x1 and x2
   if ( xL( Mode ) == xMin )
      x2 = xR( Mode );
      %x1 = max( max( xL( Mode + 1 : M ) ), xMin );
      x1 = max( xL( Mode + 1 : M ) );

      % Begin bisection (allowing no more than MaxBIS bisections per mode)
      for J = 1 : MaxBisections
         x = x1 + ( x2 - x1 ) / 2;
         [ ~, ModeCount ] = pekeris( x );
         NZeros = ModeCount - NZero1;
         
         if ( NZeros < Mode )   % not too many zeros, this is a new right bdry
            x2 = x;
            xR( Mode ) = x;
         else                        % this is a new left bdry
            x1 = x;
            if ( xR( NZeros + 1 ) >= x )
               xR( NZeros + 1 ) = x;
            end
            
            if ( xL( NZeros     ) <= x )
               xL( NZeros     ) = x;
            end
         end
         
         % when we have replaced the default, initial values, we are done
         if ( xL( Mode ) ~= xMin )
            break   % next Mode
         end
      end
   end
end % next Mode

return
