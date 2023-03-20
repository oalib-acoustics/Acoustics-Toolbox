function Slice = CreateSliceStructure( )

% Create a new Slice struct with default field values.
%  It has no input arguments and returns a Slice struct with all required
%  fields populated with default values. (A notional BELLHOP calculation
%  with Pekeris sound speed profile and water depth of 1000 meters).

% Construct the various sub-structs first

Top.HS  = struct( [] );
Top.Opt = 'CVFT  ';	% 'C' C-linear interpolation of sound speed profile
			% 'V' Vacuum above top, air-sea interface (open ocean)
			% 'F' Frequency dependent attenuation, dB/(kmHz)
			% 'T' Thorp volume attenuation formula
			% ' ' Flat air-sea interface (no *.ati altimetry file)
			% ' ' Trace all the rays in the beam fan

HS.alphaR = 1600.0;
HS.betaR  = 0.00;
HS.rho    = 1.50;
HS.alphaI = 0.75;
HS.betaI  = 0.00;

Bot.alphaR = 1600.0;
Bot.betaR  = 0.00;
Bot.rho    = 1.50;
Bot.alphaI = 0.75;
Bot.betaI  = 0.00;
Bot.HS     = HS;
Bot.Opt    = 'A*';	% Acousto-Elastic bottom, read in *.bty file

Bdry.Top = Top;
Bdry.Bot = Bot;

Box.r =   10.1;			% Maximum range for tracing rays (km)
Box.z = 1001.0;			% Maximum depth for tracing rays  (m)

Beam.alpha   = [ -45.0, 45.0 ];	% Beam fan angles (angles < 0, toward surface)
Beam.Box     = Box;
Beam.deltas  = 0.0;
Beam.Nbeams  = 0;		% Number of rays in the source beam fan
Beam.RunType = 'CB R';		% Coherent TL calculation, Gaussian beams

cInt.Low  = 1500.0;
cInt.High = 1.0d+09;

Plane.specify_method = 'POLAR';
Plane.min_r_km       = 0.0;      % Minimum range of receivers in km
Plane.max_r_km       = 10.0;     % Maximum range of receivers in km
Plane.min_d_m        = 0.0;      % Minimum depth of receivers in meters
Plane.max_d_m        = 1000.0;	% Maximum depth of receivers in meters
Plane.nrd            = 50;       % Number of receivers in depth dimension
Plane.nrr            = 40;       % Number of receivers in range dimension

r.depth = Plane.max_d_m  * ( 1 : Plane.nrd ) / Plane.nrd;
r.range = Plane.max_r_km * ( 1 : Plane.nrr ) / Plane.nrr;
s.depth = 50.0;

Pos.Nrd = Plane.nrd;
Pos.Nrr = Plane.nrr;
Pos.Nsd = 1;
Pos.r   = r;
Pos.s   = s;

RMax = 10.1;

raw( 1 ).z      = [    0.0 1000.0 ];	% Depth, values > 0 are underwater
raw( 1 ).alphaR = [ 1500.0 1500.0 ];	% P-wave phase velocity (m/s)
raw( 1 ).betaR  = [    0.0    0.0 ];	% S-wave phase velocity (m/s)
raw( 1 ).rho    = [    1.0    1.0 ];	% Fluid density (g/cm^3)
raw( 1 ).alphaI = [    0.0    0.0 ];	% P-wave attenuation (Top.Opt(3) units)
raw( 1 ).betaI  = [    0.0    0.0 ];	% S-wave attenuation (Top.Opt(3) units)

SSP.NMedia = [ 1 ];
SSP.N      = [ 0 ];
SSP.sigma  = [ 0.0 ];
SSP.depth  = [ 0.0 1000.0 ];
SSP.raw    = raw;

% Now assemble everything into a Slice struct

Slice.acoustic_model = 'BELLHOP';
Slice.bathy_delta_km = 0.05;
Slice.bathy_interp   = 'C';
Slice.bathy_rd       = [];
Slice.Bdry           = Bdry;
Slice.Beam           = Beam;
Slice.cInt           = cInt;
Slice.freq           = 5000.0;
Slice.Plane          = Plane;
Slice.Pos            = Pos;
Slice.RMax           = RMax;
Slice.SSP            = SSP;
Slice.SSP_mat        = [];
Slice.SSP_range      = [];
Slice.TitleEnv       = 'Default Slice';
