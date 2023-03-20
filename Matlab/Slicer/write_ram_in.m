function write_ram_in( TitleEnv, freq, bathy_rd, SSP, Bdry, Pos )

% WRITE_RAM_IN - write the 'ram.in' input file for use by NRL RAM

% mbp 2009
% note: I had an unusual case with a parabolic mirror where round-off in
% the receiver ranges was a problem (+/- 1 m in range). Changed the output
% format from %6.2f to %6f to accommodate that. Similar issues may occur
% elsewhere in the code below ...
%
% Edited by various other people
% variable name conventions are not consistent with use elsewhere in the
% Acoustics Toolbox (bathy); fix later (mbp)
% Also, need to check the manner in which step sizes are adjusted.

% constants

nseg_lambda = 15;	% number of segments per wavelength for depth spacing
dr_factor   = 4;	% depth spacing multiplier for range spacing
fuzz        = 1.125;	% wiggle room for output spacing

% open the file for writing

ram_in = 'ram.in';	% RAM is "hardwired" to look for this file

fid = fopen( ram_in, 'wt' );

if ( fid < 0 )
  error( [ mfilename, ': unable to open ram.in for writing' ] );
end

% derive the input values that RAM wants

zs = Pos.s.z( 1 );	% depth of the source 
zr = Pos.r.z( 1 );	% depth of the TL line

% some conservative values for the required step size ...

lambda = 1500.0 / freq;
est_dz = lambda / nseg_lambda;	% lambda / number of segments per wavelength
est_dr = dr_factor * est_dz;	% related to the amount of variation in range

dz = est_dz;
dr = est_dr;

% compute range skip factor

d_range = 1000.0 * Pos.r.r( 1 );	% convert km -> m
ndr = round( d_range / dr );
ndr = max( ndr, 1 );

% adjust dr to get consistent spacing of the output
dr = d_range / ndr;

% sanity check on the adjusted value
if ( dr > fuzz * est_dr )
  warning( [ mfilename, ': adjusted value of dr could be too coarse' ] );
end

% compute depth skip factor

d_depth = Pos.r.z( 1 );
ndz = round( d_depth / dz );
ndz = max( ndz, 1 );

% adjust dz to get consistent spacing of the output
dz = d_depth / ndz;

% sanity check on the adjusted value
if ( dz > fuzz * est_dz )
  warning( [ mfilename, ': adjusted value of dz could be too coarse' ] );
end

zmax = 1.5 * SSP.depth( 2 );			% maximum depth (m)
rmax = 1000.0 * Pos.r.r( end );		% maximum range (m)
rmax = rmax + ndr * dr;	% kludge: RAM does not always output the last range

zmplt = Pos.r.z( end );

c0 = min( SSP.raw( 1 ).alphaR );	% reference sound speed

np = 8;		% number of Pade terms
ns = 1;		% number of stability constraints
rs = 0.0;	% maximum range of stability constraints

% write the header

fprintf( fid, '''%s''\n', TitleEnv );
fprintf( fid, '%8.2f %8.2f %8.2f\t\t! freq(Hz) zs zr\n', freq, zs, zr );
fprintf( fid, '%8.2f %10.6f %4i\t\t! rmax dr ndr\n', rmax, dr, ndr );
fprintf( fid, '%8.2f %10.6f %4i %8.2f\t! zmax dz ndz zmplt\n', zmax, dz, ndz, zmplt );
fprintf( fid, '%8.2f %2i %2i %8.2f\t\t\t! c0 np ns rs\n', c0, np, ns, rs );

% write bathymetry

if ~isempty( bathy_rd )
  nBtypts = size( bathy_rd, 1 );
  bathy_r = 1000.0 * bathy_rd( :, 1 );		% convert from km -> m
  bathy_d =          bathy_rd( :, 2 );
else
  % defaults
  nBtypts = 2;
  bathy_r = [            0.0;          1.0e9 ];
  bathy_d = [ Bdry.Bot.depth; Bdry.Bot.depth ];
end

for ii = 1 : nBtypts - 1
  fprintf( fid, '%8.2f %8.2f\t\t\t! rb zb\n', bathy_r( ii ), bathy_d( ii ) );
end

fprintf( fid, '-1 -1\n' );   % terminator

% profile updates

for ii = 1 : nBtypts - 1

  if ( ii > 1 )
    fprintf( fid, '\t%8.2f\t! rp\n', bathy_r( ii ) );   % range of profile update
  end

  % sound speed in water column

  medium = 1;
  for jj = 1 : length( SSP.raw( medium ).z )
    fprintf( fid, '\t%8.2f %8.2f\t! z c\n', ...
                 SSP.raw( medium ).z( jj ), SSP.raw( medium ).alphaR( jj ) );
  end
  fprintf( fid, '-1 -1\n');   % terminator

  % sediment properties:
  
  % sound speed

  fprintf( fid, '\t%8.2f %8.2f /\t! lower halfspace\n', 0.0, Bdry.Bot.HS.alphaR );
  fprintf( fid, '-1 -1\n');   % terminator

  % density

  fprintf( fid, '\t%8.2f %8.2f\t! lower halfspace\n', 0.0, Bdry.Bot.rho );
  fprintf( fid, '-1 -1\n');   % terminator

  % attenuation

  fprintf( fid, '\t%8.2f %8.2f\t! lower halfspace\n', 0.4 * SSP.depth( 2 ), Bdry.Bot.HS.alphaI );
  fprintf( fid, '\t%8.2f %8.2f\t! lower halfspace\n', 0.5 * SSP.depth( 2 ),                 10 );
  fprintf( fid, '-1 -1\n');   % terminator

end

fclose( fid );
