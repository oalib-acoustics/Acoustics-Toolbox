
% run kraken for every environmental file

runkraken = which( 'kraken.exe' );

%%
DirInfo = dir( '*.env' );

for ifile = 1 : length( DirInfo )
    [ ~, filename, ~ ] = fileparts( DirInfo( ifile ).name );
    disp( filename )
    eval( [ '! "' runkraken '" ' filename ] );
end

%%
depths = 36:2:52;
for depth = depths
    filename = [ 'lanta' num2str( depth ) ];
    disp( filename )
    eval( [ '! "' runkraken '" ' filename ] );
end

depths = 36:2:48;
for depth = depths
    filename = [ 'lantb' num2str( depth ) ];
    disp( filename )
    eval( [ '! "' runkraken '" ' filename ] );
    
    filename = [ 'lantc' num2str( depth ) ];
    disp( filename )
    eval( [ '! "' runkraken '" ' filename ] );
    
    filename = [ 'lantd' num2str( depth ) ];
    disp( filename )
    eval( [ '! "' runkraken '" ' filename ] );
end

depths = 2:2:52;
for depth = depths
    filename = [ 'lante' num2str( depth, '%02i' ) ];
    disp( filename )
    eval( [ '! "' runkraken '" ' filename ] );
end
