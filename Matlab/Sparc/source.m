function source( T, S, ~, Nsd, Ntout, omega, ~, ~, Pulse, PulseTitle, IniFlag )

% !!!! looks like this routine is not currently used by sparcM

% Evaluate the source time series

% stemp is a temporary workspace

% IniFlag, used to control filtering, has to be controlled on the outside
% since SPARC changes the filtering as a function of wavenumber.

MaxNt  = 1000000;
MaxNST = 1000000;

%SAVE TF, SF, Nt

% If first call, tabulate the time series

if ( IniFlag )
   
   if ( ~exist( SF, 'var' ) )
      ALLOCATE( SF( Nsd, MaxNST / Nsd ) )
      if ( Pulse(1:1) == 'F' || Pulse(1:1) == 'B' )
         % From a file
         %CALL STSHDR( PulseTitle, SD, Nsd )
         %CALL SFILE( Pulse, PulseTitle, Nsd, TF, stempR, SF, Nt, MaxNt, MaxNST )
         deltat = TF( 2 ) - TF( 1 );
      else
         % Use one of the canned wavelets
         Nsd    = 1;
         deltat = 0.2 / omega;
         TSTART = -200.0 * deltat;
         Nt     = 1024;
         for it = 1 : Nt
            TF( it ) = TSTART + ( it - 1 ) * deltat;
            cans( TF( it ), omega, Pulse, SF, Nsd, it, PulseTitle )
         end
      end
      
      disp( 'Warning: SparcM does not implement the source time series filtering' )
      % Filter the time series
      %FILTER( Pulse, deltat, stemp, SF, Nsd, Nt, fmin, fmax )
      
      IniFlag = 0;
   end
   
   evalu( T, Pulse, TF, SF, Nsd, Nt, S, Ntout );
end