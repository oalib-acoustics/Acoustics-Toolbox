

% run kraken for every environmental file

depths = 36:2:52;

for ii = 1 : length( depths )
   fid = fopen( [ 'lanta' num2str( depths( ii ) ) '.env' ] );
   title = fgetl( fid );
   title = fgetl( fid );
   title = fgetl( fid );
   title = fgetl( fid );
   title = fgetl( fid );
   title = fgetl( fid );
   c( 1, depths( ii ) ) = str2num( title( 9:16 ) );  
end

depths = 36:2:48;
for ii = 1 : length( depths )
   fid = fopen( [ 'lantb' num2str( depths( ii ) ) '.env' ] );
   title = fgetl( fid );
   title = fgetl( fid );
   title = fgetl( fid );
   title = fgetl( fid );
   title = fgetl( fid );
   title = fgetl( fid );
   c( 2, depths( ii ) ) = str2num( title( 9:16 ) );
   
   fid = fopen( [ 'lantc' num2str( depths( ii ) ) '.env' ] );
   title = fgetl( fid );
   title = fgetl( fid );
   title = fgetl( fid );
   title = fgetl( fid );
   title = fgetl( fid );
   title = fgetl( fid );
   c( 3, depths( ii ) ) = str2num( title( 9:16 ) );
   
   fid = fopen( [ 'lantd' num2str( depths( ii ) ) '.env' ] );
   title = fgetl( fid );
   title = fgetl( fid );
   title = fgetl( fid );
   title = fgetl( fid );
   title = fgetl( fid );
   title = fgetl( fid );
   c( 4, depths( ii ) ) = str2num( title( 9:16 ) );
end


depths = 2:2:52;
for ii = 1 : length( depths )
   fid = fopen( [ 'lante' num2str( depths( ii ), '%02i' ) '.env' ] );
   title = fgetl( fid );
   title = fgetl( fid );
   title = fgetl( fid );
   title = fgetl( fid );
   title = fgetl( fid );
   title = fgetl( fid );
   c( 5, depths( ii ) ) = str2num( title( 9:16 ) );
end
