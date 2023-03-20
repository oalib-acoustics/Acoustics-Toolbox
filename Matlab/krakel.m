function krakel( filename )

% run the KRAKEL program
%
% usage: krakel( filename )
% where filename is the environmental file (without the extension)

runkrakel = which( 'krakel.exe' );

if ( isempty( runkrakel ) )
   error( 'krakel.exe not found in your Matlab path' )
else
   eval( [ '! "' runkrakel '" ' filename ] );
end

runfield = which( 'field.exe' );
eval( [ '! "' runfield '" ' filename ' < field.flp > field.prt' ] );

