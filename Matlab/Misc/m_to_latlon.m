function [ lat, lon ] = m_to_latlon( x, y, lat_O, lon_O )
%
% usage [ lat, lon ] = m_to_latlon( x, y, lat_O, lon_O )
% converts an (x, y) coordinate in meters
% to latitude/longitude coordinate given in decimal degrees
% relative to a coordinate origin at ( lat_O, lon_O )
% lat, lon can be scalar, vectors, or matrices
%
% mbp Feb. 2003

R_earth = 6371010.0  ; % mean radius of earth in meters (accurate to +/- 20 m
deg2rad = pi / 180;   % multiplier to convert degrees to radians

lat_O2 = deg2rad * lat_O;
lon_O2 = deg2rad * lon_O;

lat2 = y / R_earth + lat_O2;
lon2 = x ./ ( R_earth .* cos( lat2 ) ) + lon_O2;

lat   = lat2 / deg2rad;
lon   = lon2 / deg2rad;
