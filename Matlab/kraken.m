function kraken( filename )

% run the KRAKEN program
%
% usage: kraken( filename )
% where filename is the environmental file (without the extension)

runkraken = which( 'kraken.exe' );
 
if ( isempty( runkraken ) )
   error( 'kraken.exe not found in your Matlab path' )
else
   eval( [ '! "' runkraken '" ' filename ] );
end

% run field to sum the modes and compute the pressure field

% Fortran version
%runfield = which( 'field.exe' );
%eval( [ '! "' runfield '" ' filename ] );

% Matlab version
field( filename );
