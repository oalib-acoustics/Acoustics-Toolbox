function make_bdry( interp_type )

% generate Elliptic bottom bathymetry
% Based on Diana McCammon's test case and her code
% Generates r,z vectors with a constant chord length
% mbp Feb. 2009

a = 20177.67;
b = 10044.32;

h = 17500.;
r1 = 0.;
z1 = 5000.;

gam = 40.;   % chord length

rmax = 37600.;

np = 1000;
eps = 0.005;
eps = 0.0000001;

z = np;
r = np;
chord = np;   % mbp: not clear that this does anything
r( 1 ) = r1;
z( 1 ) = z1;

for i = 2 : np
   % finding each point
   rn    = r ( i - 1 ) + gam;   % initial guess to rn
   delta = 100.;
   
   % Newton iteration for solution to combinational equation
   
   while abs( delta ) > eps
      fofr  = ( b * sqrt( 1 - ( rn - h )^2 / a^2 ) - z( i - 1 ) )^2 + ( rn - r( i - 1 ) )^2 - gam^2;
      fprim = b / 2 * 1. / sqrt( 1. - ( rn - h )^2 / a^2 ) + rn;
      delta = -fofr / fprim;
      rn    = rn + delta;
   end
   
   r( i )         = rn;
   z( i )         = b * sqrt( 1. - ( r( i )- h )^2 / a^2 );
   chord( i - 1 ) = sqrt( ( r( i ) - r( i - 1 ) ) ^2 + ( z( i ) - z( i - 1 ) )^2 );
   
   % jump out if we get a range outside rmax
   if r( i ) > rmax
      break
   end
end

% write the alitmetry file
atifil = 'Ellipse.ati';
rngdep = [ r'/1000 -z' ];
writebdry( atifil, interp_type, rngdep )

% write the bathymetry file
btyfil = 'Ellipse.bty';
rngdep = [ r'/1000 z' ];
writebdry( btyfil, interp_type, rngdep )
