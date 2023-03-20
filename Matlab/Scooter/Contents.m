% SCOOTER
%
% Files
%   backsub     - Performs back-substitution for a symmetric tridiagonal linear system
%   backsubOld  - x = backsub( N, mults, d, e, b )
%   bcimp       - Compute Boundary Condition IMPedance
%   bounceM     - Wavenumber integration code for ocean acoustics problems
%   elasdn      - Propagates through an elastic layer using
%   elasup      - Propagates through an elastic layer using
%   factortri   - [ mults, dt, et ] = factor( N, d, e )
%   fieldsco    - Calculates the field using the Green's function produced by SCOOTER
%   fieldscoOrg - calculates the field using the Green's function produced by SCOOTER
%   FTS         - Fourier transform of SCOOTER Green's function to produce pressure
%   FTSOrg      - Fourier transform of SCOOTER Green's function to produce pressure
%   HTS         - Hankel transform of SCOOTER Green's function to produce pressure
%   HTSOrg      - Hankel transform of SCOOTER Green's function to produce pressure
%   init_matrix - Initializes arrays defining difference equations
%   kernel      - Solve system for a sequence of k-values
%   profil_old  - return the SSP values at user specified points

%   scooterM    - Wavenumber integration code for ocean acoustics problems
%   Solve       - Set up the linear system and solve
%   weight      - Computes weights for linear interpolation
%   fieldscoFFT - Calculates the field using the Green's function produced by SCOOTER
%   backsubV    - Performs back-substitution for a symmetric tridiagonal linear system
%   factortriV  - [ mults, dt, et ] = factor( N, d, e )
