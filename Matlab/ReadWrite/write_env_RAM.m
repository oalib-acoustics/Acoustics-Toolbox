function write_env_RAM( ~, ~, TitleEnv, freq, SSP, Bdry, Pos, ~, ~, ~, varargin )

% Write an environmental file
% mbp 2009
% note: I had an unusual case with a parabolic mirror where round-off in
% the receiver ranges was a problem (+/- 1 m in range)
% Changed the output format from %6.2f to %6f to accommodate that.
% Similar issues may occur elsewhere in the code below ...

global xBot NbtyPts

envfil = 'ram.in';
btyfil = 'temp.bty';
BotBTY = '*';
depthB = 1e9;
rBox   = 1e9;
RMax   = Pos.r.r( end );

fid = fopen( envfil, 'wt' );   % create new envfil

zmax = 1.5 * SSP.depth( 2 );   % extend water depth by 50%

% some conservative guesses at the required step size ...
dz = 1500 / freq / 15;   % lamda / 15
dr = 10 * dz;            % this one should be related to the amount of variation in range

% dr = 1000 * ( Pos.r.r( end ) - Pos.r.r( 1 ) ) / length( Pos.r.r );   % range step
% dz =        ( Pos.r.z( end ) - Pos.r.z( 1 ) ) / length( Pos.r.r );   % depth step

% round the mesh spacing to be an even divider of the water depth
% nn = round( zmax / dz );
% dz = zmax / nn;
% nn = round( 1001 * RMax / dr ) - 1;
% dr = 1001 * RMax / nn;
% 
% ndr = 1;   % range skip factor
% ndz = 1;   % depth skip factor

ndr = round( ( Pos.r.r( 2 ) - Pos.r.r( 1 ) ) / dr );   % range skip factor
ndz = round( ( Pos.r.z( 2 ) - Pos.r.z( 1 ) ) / dz );   % depth skip factor

ndr = max( ndr, 1 );
ndz = max( ndz, 1 );

zmplt = Pos.r.z( end );

c0 = mean( SSP.raw( 1 ).alphaR );   % reference sound speed

np = 4;   % number of Pade terms
ns = 4;   % number of starter terms
rs = 0.0; % range of source?

% header

fprintf( fid, '''%s'' ! Title', TitleEnv );
fprintf( fid, '%8.2f %8.2f %8.2f \t \t \t ! Frequency (Hz) zs zr', freq, Pos.s.z( 1 ), Pos.r.z( end ) );
fprintf( fid, '%8.2f %8.2f %3i   \t \t \t ! rmax dr ndr', 1000 * RMax, dr, ndr );
fprintf( fid, '%8.2f %8.2f %3i %8.2f   \t \t \t ! zmax dz ndz zmplt', zmax, dz, ndz, zmplt );
fprintf( fid, '%8.2f %1i %1i %8.2f   \t \t \t ! c0 np ns rs', c0, np, ns, rs );

% bathymetry
BotBTY = Bdry.Bot.Opt( 2 : 2 );

if ( exist( btyfil, 'file' ) && BotBTY == '*' )
   readbty( btyfil, BotBTY, depthB, rBox )
else
   NbtyPts = 3;
   clear xBot
   xBot( 1, : ) = [ 0 0 1e9 ];
   xBot( 2, : ) = [ Bdry.Bot.depth Bdry.Bot.depth Bdry.Bot.depth ];
end

for ii = 2 : NbtyPts - 1
   fprintf( fid, '%8.2f %8.2f   \t \t \t ! rb zb', xBot( 1, ii ), xBot( 2, ii ) );
end
fprintf( fid, '-1 -1' );   % terminator

% SSP
for ii = 2 : NbtyPts - 1
   
   if ( ii ~= 2 )
      fprintf( fid, '\t %6.2f \t ! rp', xBot( 1, ii ) );   % range for this profile
   end
   
   % SSP
   
   medium = 1;
   for jj = 1 : length( SSP.raw( medium ).z )
      fprintf( fid, '\t %6.2f %6.2f \t ! z c', ...
         [ SSP.raw( medium ).z( jj ) SSP.raw( medium ).alphaR( jj ) ] );
   end
   fprintf( fid, '-1 -1' );   % terminator
   
   % Sediment
   fprintf( fid, '    %6.2f %6.2f /  \t ! lower halfspace', 0.0, Bdry.Bot.HS.alphaR );
   fprintf( fid, '-1 -1' );   % terminator
   
   % density
   fprintf( fid, '    %6.2f %6.2f  \t ! lower halfspace', 0.0, Bdry.Bot.rho );
   fprintf( fid, '-1 -1' );   % terminator
   
   % attenuation
   fprintf( fid, '    %6.2f %6.2f  \t ! lower halfspace', 0.0, Bdry.Bot.HS.alphaI );
   fprintf( fid, '    %6.2f %6.2f  \t ! lower halfspace', 1.9 * SSP.depth( SSP.NMedia+1 ), Bdry.Bot.HS.alphaI );
   fprintf( fid, '    %6.2f %6.2f  \t ! lower halfspace', 2.0 * SSP.depth( SSP.NMedia+1 ), 10 );
   fprintf( fid, '-1 -1' );   % terminator

end

fclose( fid );

