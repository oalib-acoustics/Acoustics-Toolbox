function [ cp, cs,rho ] = profil( ~, Med, N1, ~, ~, ~ )

% return the SSP values at user specified points

ILoc = Loc( Med );
N    = N1 - 1;

ZT   = linspace( Z( ILoc ), Z( ILoc + 1 ), N );
cp   = interp( ZT, cpv  );
cs   = interp( ZT, csv  );
rho  = interp( ZT, rhov );

