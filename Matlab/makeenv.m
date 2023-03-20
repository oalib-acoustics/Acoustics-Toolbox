% makeenv
% set up the workspace variables with defaut values for an envfil 
% michael b. porter October 23, 1999

clear all

envfil = 'test.env';

title  = '''sample''';
freq   = 50.0;
Nmedia = 1;

topopt = '''NVF''';

Npts  = 0;
sigma = 0;
Depth = 5000.0;

% sound speed profile

cofz = [ 0 1500;
  5000 1500 ];

% lower halfspace

botopt = '''A''';

cpb  = 2000;
csb  = 0;
rhob = 2.0;

cmin = 1400;
cmax = 2000;

% source/receiver depths

sd  = 500;
Nsd = length( sd );

rd  = 2500;
Nrd = length( rd );

Nrr = 2;
rr  = [ 100.0 200.0 ];

runtyp = '''R''';

Nbeams = 80;
alpha  = [ -20 20 ];

step = 0.0;
zbox = 1.1 * Depth;
rbox = 10.0;
