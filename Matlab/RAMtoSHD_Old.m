function RAMtoSHD( SHDFile )
%
% RAMtoSHD
% Converts the standard RAM output file to a shdfil
% the SHDFile must be a complete file name, including the extension

filename = 'tl.grid';
fid = fopen( filename, 'rb' );

if ( fid == -1 )
   error( 'No shade file with that name exists; you must run a model first' );
end

%
% read the header with the problem description
%   write( 3 ) freq,zs,zr,rmax,dr,ndr,zmax,dz,ndz,zmplt,c0,np,ns,rs
%
% FORTRAN standards do not specify precisely what gets
% written to a file when using unformatted write statements (so this
% is compiler-dependent).
%
% Often times, the record starts with an integer word that gives the
% data length (sometimes in bytes, sometimes in words). The code below
% is for gfortran, which starts AND ends the record with an integer
% word with the data length in units of bytes.

% read in the grid size information from the header
gfortran_lrecl    = fread( fid, 1, 'float32' );   % Most Fortran write lrecl first
dz    = fread( fid, 1, 'float32' );
ndz   = fread( fid, 1, 'int32'   );
zmplt = fread( fid, 1, 'float32' );
dr    = fread( fid, 1, 'float32' );
ndr   = fread( fid, 1, 'int32'   );
rmax  = fread( fid, 1, 'float32' );
freq  = fread( fid, 1, 'float32' );
zs    = fread( fid, 1, 'float32' );
fseek( fid, 56, 'bof' ); %reposition to end of first record

% Pos vectors
% Note that the range in RAM may not write a value at rmax
% if rmax is not an integer multiple of ndr * dr
% In that case the vector stops at an earlier range
Pos.s.depth = zs;
Pos.r.depth = 0  : ndz * dz : zmplt; % + ndz * dz;
Pos.r.range = dr : ndr * dr : rmax;

Nrd = length( Pos.r.depth );
Nrr = length( Pos.r.range );

TL = fread( fid, [ Nrd, Nrr ], 'float32' );    %Read complex data

fclose( fid );

PlotTitle = 'RAM';
PlotType  = 'rectilin  ';
atten     = 0;
pressure( 1, 1, :, : ) = 10 .^ ( -TL / 20 );
save( SHDFile, 'PlotTitle', 'PlotType', 'freq', 'atten', 'Pos', 'pressure' )  % problem: sd is outside the loop
