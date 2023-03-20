function krakenc_nofield( filename )

% run the KRAKENC program
%
% usage: krakenc( filename )
% where filename is the environmental file


runkrakenc = which( 'krakenc.exe' );

if ( isempty( runkrakenc ) )
   error( 'krakenc.exe not found in your Matlab path' )
else
   eval( [ '! "' runkrakenc '" ' filename ] );
end

