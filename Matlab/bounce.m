function bounce( filename )

% run the BOUNCE program
%
% usage: bounce( filename )
% where filename is the environmental file

runbounce = which( 'bounce.exe' );

if ( isempty( runbounce ) )
   error( 'bounce.exe not found in your Matlab path' )
else
   eval( [ '! "' runbounce '" ' filename ] );
end

