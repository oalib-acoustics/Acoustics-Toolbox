%Calculates the distance between two geographical points.
%
%SYNOPSIS:  r     = distance( gci, gcf )
%
%          where gci and gcf corresponds to the initial and final 
%          geographical coordinates of the corresponding points. 
%          These coordinates should be given in decimal format [lat long].
%          "r" is the output value in kilometers.
%          mbp: I believe r is in meters since R_earth is in meters
%
% Either one of both of gci or gcf can be vectors
%
%See also distlat and distlong.

function  r = distance( gc1, gc2 )

% make sure lat/long are first/second columns

if ( size( gc1, 2 ) ~= 2 )
    gc1 = gc1';
end

if ( size( gc2, 2 ) ~= 2 )
    gc2 = gc2';
end

lai = gc1( :, 1 );   % initial coordinates
loi = gc1( :, 2 );

laf = gc2( :, 1 );   % final   coordinates
lof = gc2( :, 2 );
                           
Cearth = 40000000;	% circumference of earth
disty = ( laf - lai ) * Cearth / 360.0;
distx = ( lof - loi ) * Cearth / 360.0 .* cos( 0.5 * ( laf + lai ) * pi / 180 );

r = sqrt( distx .* distx + disty .* disty );
