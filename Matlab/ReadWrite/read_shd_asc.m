function [ PlotTitle, PlotType, freqVec, freq0, atten, Pos, pressure ] = read_shd_asc( filename )

% Read an ascii shade file

% open the file
fid = fopen( filename, 'r' );
if ( fid == -1 )
    error( 'No shade file with that name exists; you must run a model first' );
end

% read

PlotTitle = fgetl( fid );
PlotType  = fgetl( fid );
Nfreq     = fscanf( fid, '%i', 1 );
Ntheta    = fscanf( fid, '%i', 1 );
Nsd       = fscanf( fid, '%i', 1 );
Nrd       = fscanf( fid, '%i', 1 );
Nrr       = fscanf( fid, '%i', 1 );
freq0     = fscanf( fid, '%f', 1 );
atten     = fscanf( fid, '%f', 1 );

freqVec     = fscanf( fid, '%f', Nfreq );
Pos.theta   = fscanf( fid, '%f', Ntheta );
Pos.s.z = fscanf( fid, '%f', Nsd );
Pos.r.z = fscanf( fid, '%f', Nrd );
Pos.r.r = fscanf( fid, '%f', Nrr );

isd = 1;
for ii = 1:isd
   temp1   = fscanf( fid, '%f', [ 2 * Nrr, Nrd ] );
end

fclose( fid );

% joint real and imaginary parts into a complex matrix

pressure = temp1( 1 : 2 : 2 * Nrr, : )' + 1i * temp1( 2 : 2 : 2 * Nrr, : )';
