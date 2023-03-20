function sparc( filename )

% run the SPARC program
%
% usage: sparc( filename )
% where filename is the environmental file

runsparc = which( 'sparc.exe' );

if ( isempty( runsparc ) )
   error( 'sparc.exe not found in your Matlab path' )
else
   eval( [ '! "' runsparc '" ' filename ] );
end