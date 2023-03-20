function make_bdry( interp_type )

% generate parabolic bottom bathymetry
% mbp Dec. 27, 2005
% example from:
% R. W. McGirr, D.B. King, J.A. Davis, J. Campbell, "An evaluation of
% range-dependent ray theory models", NORDA report 115.

b = 250000;
c = 250;

r = 0 : 100 : 24750;
z = 0.002 * b * sqrt( 1 + r / c );

clear r z
%z = 0 : 25 : 5000;   % for ray plot in JKPS
z = 0 :  5 : 5000;   % for TL comparison to virtual source method
r = c * ( ( z / 0.002 / b ).^2 - 1 ) / 1000;   % range in km

%%
rngdep = [ r' -z' ];
writebdry( 'ParaBot.ati', interp_type, rngdep )

rngdep = [ r' +z' ];
writebdry( 'ParaBot.bty', interp_type, rngdep )