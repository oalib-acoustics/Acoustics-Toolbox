function Pos = readszrz( fid )

% readszrz Read source depths and receiver depths
%
% Variable 'Pos' is a structure:
%
% Pos.r.z = vector of receiver depths
% Pos.Nrz = number of receiver depths
% Pos.s.z = vector of source   depths
% Pos.Nsz = number of source   depths

%%
% source depths
fprintf( '\n_______________________ \n' )

[ Pos.s.z, Pos.Nsz ] = readvector( fid );

fprintf( '\n Number of source   depths, NSz   = %i \n', Pos.Nsz )
fprintf( '\n Source depths, Sz (m)\n' );

if ( Pos.Nsz < 10 )
   fprintf( '%8.2f \n', Pos.s.z )   % print all the depths
else
   fprintf( '%8.2f ... %8.2f \n', Pos.s.z( 1 ), Pos.s.z( end ) ) % print first, last depth
end

%%
% receiver depths
fprintf( '\n_______________________ \n' )

[ Pos.r.z, Pos.Nrz ] = readvector( fid );

fprintf( '\n Number of receiver depths, NRz   = %i \n', Pos.Nrz )
fprintf( '\n Receiver depths, Rz (m)\n' );

if ( Pos.Nrz < 10 )
   fprintf( '%8.2f \n', Pos.r.z )   % print all the depths
else
   fprintf( '%8.2f ... %8.2f \n', Pos.r.z( 1 ), Pos.r.z( end ) ) % print first, last depth
end

