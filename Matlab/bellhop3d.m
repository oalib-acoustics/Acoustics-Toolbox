function bellhop3d( filename )

% run the BELLHOP3D program
%
% usage: bellhop( filename )
% where filename is the environmental file

runbellhop3d = which( 'bellhop3d.exe' );

if ( isempty( runbellhop3d ) )
   error( 'bellhop3d.exe not found in your Matlab path' )
else
   eval( [ '! "' runbellhop3d '" ' filename ] );
end