function [    x2, Tray2, p2, q2, tau2, Rfa2, c2 ] = ...
   step( SSP, x0, Tray0, p0, q0, tau0, Rfa0, ...
   xTop, nTop, xBot, nBot,  rTopSeg, rBotSeg, deltas )

% Does a single step along the ray
% x denotes the ray coordinate, (r,z)
% Tray denotes the scaled tangent to the ray

global Layer

% *** Phase 1 of modified polygon method (an Euler step) ***

[ c0, c0_imag, gradc0, crr0, crz0, czz0, Layer ] = ssp( x0, SSP, Layer );

csq0 = c0 * c0;
cnn0_csq0 = crr0 * Tray0(2)^2 - 2.0 * crz0 * Tray0(1) * Tray0(2) + czz0 * Tray0(1)^2;

Layer0 = Layer;   % make note of current layer

h = deltas;   % initially set the step h, to the basic one, deltas
h = reducestep( SSP, x0, Tray0, Layer0, c0, xTop, nTop, xBot, nBot,  rTopSeg, rBotSeg, deltas, h );
halfh = 0.5 * h;   % first step of the modified polygon method is a half step

x1    = x0    + halfh * c0 * Tray0;
Tray1 = Tray0 - halfh * gradc0 / csq0;
p1    = p0    - halfh * cnn0_csq0 * q0;
q1    = q0    + halfh * c0        * p0;

% *** Phase 2 of modified polygon method ***

[ c1, c1_imag, gradc1, crr1, crz1, czz1, Layer ] = ssp( x1, SSP, Layer );

csq1 = c1 * c1;
cnn1_csq1 = crr1 * Tray1(2)^2 - 2.0 * crz1 * Tray1(1) * Tray1(2) + czz1 * Tray1(1)^2;

h = reducestep( SSP, x0, Tray1, Layer0, c1, xTop, nTop, xBot, nBot,  rTopSeg, rBotSeg, deltas, h );

% use blend of f' based on proportion of a full step used.
w1 = h / ( 2.0d0 * halfh );
w0 = 1.0d0 - w1;

% you could absorb h into w0, w1 but I leave as follows for clarity, hoping compiler is smart enough to optimize
x2     = x0    + h * ( w0 * c0 * Tray0     + w1 * c1 * Tray1     );
Tray2  = Tray0 - h * ( w0 * gradc0 / csq0  + w1 * gradc1 / csq1  );
p2     = p0    - h * ( w0 * cnn0_csq0 * q0 + w1 * cnn1_csq1 * q1 );
q2     = q0    + h * ( w0 * c0        * p0 + w1 * c1        * p1 );

c0_cmplx = c0 + 1i * c0_imag;
c1_cmplx = c1 + 1i * c1_imag;
tau2   = tau0  + h * ( w0 / c0_cmplx       + w1 / c1_cmplx       );

Rfa2   = Rfa0;

% If we crossed an interface, apply jump condition

[ c2, ~, gradc2, ~, ~, ~, Layer ] = ssp( x2, SSP, Layer );

if ( Layer ~= Layer0 )
   RN  = -Tray2( 1 )^2 / Tray2( 2 ) * ( gradc2( 2 ) - gradc0( 2 ) ) / c0;
   % above needs updating for c(r,z) problem
   p2 = p2 + q2 * RN;
end
