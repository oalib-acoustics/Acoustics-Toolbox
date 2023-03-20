function AddArr( omega, id, ir, amp, delay, angle, ~ )

% Add the amplitude and delay for an ARRival into a matrix of same.
% Extra logic to keep only the strongest arrivals

% This code needs lots of work; under construction

global Arr

PhaseTol = 0.5;   % arrivals with essential the same phase are grouped into one

Nt     = Arr.Narr( id, ir );    % # of arrivals
NewRay = 1;

% Is this the second bracketting ray of a pair?
% If so, we want to combine the arrivals to conserve space.
% (test this by seeing if the arrival time is close to the previous one)
% (also need that the phase is about the same to make sure surface and direct paths are not joined)

% search for old arrivals which are possible bracketting pairs
if ( Nt >= 1 )
   % search for entries in arrivals table with near same arrival time
   if ( omega * abs( delay - Arr.delay( id, ir, Nt ) ) < PhaseTol && ...
         abs( Arr.phase( id, ir, Nt ) - Phase ) < PhaseTol )
      NewRay = 0;
   end
end

if ( NewRay )
   if ( Nt >= MxNarr )   % space available to add an arrival?
      [ ~, iArr ] = min( Arr.A( id, ir, : ) ); % no space available; substitute new arrival for weakest one
      if ( amp > Arr.A( id, ir, iArr ) )
         Arr.A(     id, ir, iArr ) = amp;  	 % amplitude
         Arr.delay( id, ir, iArr ) = delay;	 % delay time
         Arr.angle( id, ir, iArr ) = angle;	 % angle
      end
      
   else     % space available; add another arrival to the table
      Arr.Narr(  id, ir         ) = Nt + 1;   % # of arrivals
      Arr.A(     id, ir, Nt + 1 ) = amp;		 % amplitude
      Arr.delay( id, ir, Nt + 1 ) = delay;	 % delay time
      Arr.angle( id, ir, Nt + 1 ) = angle;	 % angle
   end
else   % second arrival of a bracketting pair; add in to existing slot
   ampTot = Arr.A( id, ir, Nt ) + amp;
   Arr.A(     id, ir, Nt ) =           Arr.A(   id, ir, Nt ) + amp;
   Arr.delay( id, ir, Nt ) = ( Arr.A( id, ir, Nt ) * Arr.delay( id, ir, Nt ) + amp * delay ) / ampTot;   % weighted sum
   Arr.ang(   id, ir, Nt ) = angle;
end
