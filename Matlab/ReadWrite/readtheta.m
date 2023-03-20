function [ theta ] = readtheta( fid )

% Read receiver angles
fprintf( '\n_______________________ \n' )
[ theta, Ntheta ] = readvector( fid );

fprintf( '\n Number of receiver angles = %i \n', Ntheta )
fprintf( '\n Receiver angles (degrees) \n' )
fprintf( '%f ', theta )
fprintf( '\n' )
