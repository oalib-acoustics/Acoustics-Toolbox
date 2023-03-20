function [ Rr ] = readr( fid )

% Read receiver ranges
fprintf( '\n_______________________ \n' )

[ Rr, NRr ] = readvector( fid );

fprintf( '\n Number of receiver ranges, NRr = %i \n', NRr )
fprintf( '\n Receiver ranges, Rr (km) \n' )

Rr = Rr';

if ( NRr < 10 )
   fprintf( '%8.2f \n', Rr )   % print all the ranges
else
   fprintf( '%8.2f ... %8.2f \n', Rr( 1 ), Rr( end ) ) % print first, last range
end
