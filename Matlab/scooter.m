
function scooter( filename )

% run the SCOOTER program
%
% usage: scooter( filename )
% where filename is the environmental file

runscooter = which( 'scooter.exe' );

if ( isempty( runscooter ) )
   error( 'scooter.exe not found in your Matlab path' )
else
   eval( [ '! "' runscooter '" ' filename ] );
end

% Fortran fields routine
% This version is more efficient in its memory useage
% It is not maintained.

%runfields = which( 'fields.exe' );
%eval( [ '! "' runfields '" ' filename ] );

% Matlab fields routine
% Much simpler code and, because of auto-parallization in Matlab is often faster
% It also allows an irregular grid of receiver ranges since it does a DFT rather than an FFT

fieldsco( [ filename '.grn' ] );
