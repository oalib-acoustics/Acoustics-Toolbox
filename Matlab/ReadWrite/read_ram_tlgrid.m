function [ PlotTitle, PlotType, freq, atten, Pos, pressure ] = read_ram_tlgrid( )

% Read the binary tl.grid output file written by RAM

% open the sequential, unformatted tl.grid file written by RAM

fid = fopen('tl.grid', 'rb', 'native');

if ( fid < 0 )
  error( [ mfilename, ': unable to open tl.grid' ] );
end

% NONE of the FORTRAN language standards (66, 77, 90, 95, 2000, 2003) specify
% the precise format of sequential, unformatted files. This makes the task
% of reading them from other languages a dicey (and sometimes frustrating)
% proposition. The code below attempts to determine the datatype used for the
% record markers. (This is the most commonly encountered difference among
% compilers, and even compiler releases!)


% first guess: unsigned 8 byte integer
rlen = fread( fid, 1, 'ubit64' );

% check the value read against the known value
if rlen == 60
  rm_datatype = 'ubit64';
else
  % try unsigned 4 byte integer
  frewind( fid );
  rlen = fread( fid, 1, 'ubit32' );
  % check the value read against the known value
  if rlen == 60
    rm_datatype = 'ubit32';
  else
    fclose( fid );
    error( [ mfilename, ': error while determining datatype of record markers' ] );
  end
end

% continue reading the header containing the problem description

freq  = fread( fid, 1, 'float' );
zs    = fread( fid, 1, 'float' );
zr    = fread( fid, 1, 'float' );
rmax  = fread( fid, 1, 'float' );
dr    = fread( fid, 1, 'float' );
ndr   = fread( fid, 1, 'int32' );
zmax  = fread( fid, 1, 'float' );
dz    = fread( fid, 1, 'float' );
ndz   = fread( fid, 1, 'int32' );
zmplt = fread( fid, 1, 'float' );
c0    = fread( fid, 1, 'float' );
np    = fread( fid, 1, 'int32' );
ns    = fread( fid, 1, 'int32' );
rs    = fread( fid, 1, 'float' );
lz    = fread( fid, 1, 'int32' );
rlen  = fread( fid, 1, rm_datatype);		% record markers

% read the TL matrix

lr = floor( rmax / ( dr * ndr ) );

TL = zeros( lz, lr );

for j = 1 : lr
  rlen = fread( fid, 1, rm_datatype );		% record markers
  TL( :, j ) = fread( fid, lz, 'float' );
  rlen = fread( fid, 1, rm_datatype );		% record markers
end

fclose( fid );   % close file

% convert to standard Acoustics Toolbox variable names

PlotTitle = 'RAM';
PlotType  = 'rectilin';
atten     = 0;

Pos.s.z     = zs;
Pos.r.z     = [ 1 : lz ] * ndz * dz;
Pos.r.r = [ 1 : lr ] * ndr * dr;

pressure( 1, 1, 1 : lz, 1 : lr ) = 10.^( -TL / 20 );




