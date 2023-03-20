function c = Munk( z )

% Munk profile 


c0   = 1500.0;
zt   = 2.0 * ( z - 1300.0 ) / 1300.0;  % used to be called x
DxDz = 2.0 / 1300.0;
c    = c0 * ( 1.0 + 0.00737 * ( zt - 1.0 + exp( -zt ) ) );

return