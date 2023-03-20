function z = sub_seamount( parameters, x, y )
% seamount model 
% x and y need to be the same size

z1 = nan( size( x ) );

r = sqrt( ( x( : ) - parameters.tip_x( 1 )).^2 + ( y( : ) - parameters.tip_y( 1 )).^2 );
r(r(:)>=parameters.radius( 1 )) = parameters.radius( 1 ); 
z1(:) = parameters.tip_z( 1 ) + parameters.slope( 1 ) * r; 

%%
z2 = nan(size(x));

r = sqrt( ( x( : ) - parameters.tip_x( 2 )).^2 + ( y( : ) - parameters.tip_y( 2 )).^2 );
r(r(:)>=parameters.radius( 2 )) = parameters.radius( 2 ); 
z2(:) = parameters.tip_z( 2 ) + parameters.slope( 2 ) * r; 

z = min(z1,z2);



return