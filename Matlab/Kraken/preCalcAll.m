function preCalcAll( FileRoot, D )

% Calculate the field for a wedge ocean
% Can then be used for an arbitrary bathymetry
% Useage:
%    precalc( FileRoot, D)
%    FileRoot is an ENVFIL
%    D is a vector of bottom depths

runkraken = which( 'krakenc.exe' );

% read in the template environment
model  = 'KRAKENC';
[ TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, fid ] = read_env( FileRoot, model );

envfil = sprintf( '%s_rd', FileRoot );

for ir = 1 : length( D )

   % substitute using the depth from the bathymetry file
   SSPtemp = SSP;   % make copy of template for temporary changes
   SSPtemp.depth( 2 ) = D( ir );

   % remove depths that exceed the bathymetry
   ii = find( SSP.raw( 1 ).z > D( ir ), 1, 'first' );
   
   if ( ~isempty( ii ) )
      SSPtemp.raw( 1 ).z = SSP.raw( 1 ).z( 1 : ii - 1 );
   
      % add a new SSP point interpolated to the bathymetry
      if ( D( ir ) > SSP.raw( 1 ).z( ii - 1 ) )   % make sure added point is strictly greater in depth than current value
         SSPtemp.raw( 1 ).z(      ii ) = D( ir );
         SSPtemp.raw( 1 ).alphaR( ii ) = interp1( SSP.raw( 1 ).z, SSP.raw( 1 ).alphaR, D( ir ) );
      end
      
      % shift next layer up
      SSPtemp.raw( 2 ).z( 1   ) = D( ir );
   end
   
   TitleEnv = [ 'Profile #' int2str( ir ) ];
   
   if ( ir == 1 )
      write_env( envfil, model, TitleEnv, freq, SSPtemp, Bdry, Pos, Beam, cInt, RMax )
   else
      write_env( envfil, model, TitleEnv, freq, SSPtemp, Bdry, Pos, Beam, cInt, RMax, 'a' )
   end
end	% next range

eval( [ '! "' runkraken '" ' envfil ] );   % run KRAKEN/FIELD
