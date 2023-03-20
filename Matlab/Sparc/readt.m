function [ t ] = readt( fid )

% Reads in output times

Nt = fscanf( fid, '%i', 1 );
fprintf( '\nNumber of output times = %i \n', Nt )
fprintf( 'Times (s) \n' )
fgetl( fid );

t = fscanf( fid, '%f', 2 );
fprintf( '%f   ', t )
fprintf( '\n' )

if Nt > 2
  t = linspace( t( 1 ), t( 2 ), Nt )'; % generate vector of receiver ranges
end

fgetl( fid );