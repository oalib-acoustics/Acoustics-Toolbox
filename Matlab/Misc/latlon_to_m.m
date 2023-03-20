function [ x, y ] = latlon_to_m( lat, lon, lat_O, lon_O )
%
% usage [ x, y ] = lat_lon_to_meters( lat, lon, lat_O, lon_O )
% converts a latitude/longitude coordinate given in decimal degrees
% to an (x,y) coordinate in meters
% relative to a coordinate origin at ( lat_O, lon_O )
% lat, lon can be scalar, vectors, or matrices
%
% !!! caution: this code is not as accurate as it should be
% mbp Feb. 2003

R_earth = 6371010.0  ; % mean radius of earth in meters (accurate to +/- 20 m
R_equatorial = 6378100.0; % accurate to 5 digits
R_polar      = 6356800.0;

deg2rad = pi / 180;   % multiplier to convert degrees to radians

lat2   = deg2rad * lat;
lon2   = deg2rad * lon;
lat_O2 = deg2rad * lat_O;
lon_O2 = deg2rad * lon_O;

% estimate radius of earth at this latitude
R = ( 1 - abs( lat_O/90 ) ) * R_equatorial + abs( lat_O/90 ) * R_polar;

x = R * ( lon2 - lon_O2 ) .* cos( lat2 );
y = R * ( lat2 - lat_O2 );

