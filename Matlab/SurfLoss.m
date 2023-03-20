function Loss = SurfLoss( theta, U, f ) 
% surface loss
% bubble losses are assumed to be the dominant effect for the total field
% (as opposed to just the coherent field)
% see the discussion in
% High-Frequency Ocean Environmental Acoustic Models Handbook

% theta = grazing angle (degrees)
% U     = wind speed (m/s)
% f     = frequency (kHz)
% SBL   = surface bubble loss (dB)

% mbp, Oct. 2002

if U > 6
    SBL = 1.26e-3 * U^1.57 * f^0.85 / sin( theta * pi / 180 );
else
    SBL = 1.26e-3 * 6^1.57 * f^0.85 / sin( theta * pi / 180 ) * exp( 1.2 * ( U - 6 ) );
end

Loss = SBL;

