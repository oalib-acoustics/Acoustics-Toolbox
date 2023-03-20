function ram

% run the RAM program
%
% usage: ram
% where filename is the environmental file (without the extension)

runram = which( 'ram.exe' );
eval( [ '! "' runram '" ' ] );
RAMtoSHD( 'RAM.shd.mat' )
