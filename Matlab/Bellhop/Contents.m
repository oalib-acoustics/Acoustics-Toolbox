% BELLHOP
%
% Files
%   AddArr               - Add the amplitude and delay for an ARRival into a matrix of same.
%   bellhopM             - Beam-tracing code for ocean acoustics problems
%   crci                 - Convert real wave speed and attenuation to a single complex wave speed

%   makeshdarr           - Make a shade file using the arrivals information

%   trace                - Traces the beam corresponding to a particular take-off angle
%   reflect              - Reflect a ray/beam off a boundary
%   scalep               - Scale the pressure field
%   step                 - Does a single step along the ray
%   topbot               - Handles top and bottom boundary conditions
%   reducestep           - Reduces the ray step size to make sure we land on interfaces and boundaries

%   ssp                  - tabulates the sound speed profile and its derivatives
%   ssp binterp       - 
%   ssp_cubic            - tabulates the sound speed profile and its derivatives

%   Munk                 - Munk profile given an an analytic function
%   Munk_interp_tests    - 
%   Munk_rd_axis         - 
%   Munkr                - Range-dependent Munk profile given as an analytic function

%   InfluenceGeoGaussian - Computes the beam influence, i.e.
%   InfluenceGeoHat      - Computes the beam influence, i.e.
