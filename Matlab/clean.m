% deletes all temporary files in the current directory

delete *.prt
delete *.shd
delete *.grn
delete *.mod
delete *.mat
delete *.arr
delete *.asv
delete *.rts
delete *.ray
delete *.irc
delete *.brc
delete *.avi

currentFolder = cd;

fileExisting  = ( exist( fullfile( currentFolder, 'tl.grid' ), 'file' ) == 2);
if ( fileExisting ); delete tl.grid; end

fileExisting  = ( exist( fullfile( currentFolder, 'tl.line' ), 'file' ) == 2);
if ( fileExisting ); delete tl.line; end

%fileExisting  = ( exist( fullfile( currentFolder, 'clean.m' ), 'file' ) == 2);
%if ( fileExisting ); clean; end

delete ._*

