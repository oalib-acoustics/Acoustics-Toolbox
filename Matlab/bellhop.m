function bellhop( filename )

% run the BELLHOP program
%
% usage: bellhop( filename )
% where filename is the environmental file

runbellhop = which( 'bellhop.exe' );

if ( isempty( runbellhop ) )
   error( 'bellhop.exe not found in your Matlab path' )
else
   eval( [ '! "' runbellhop '" ' filename ] );
end