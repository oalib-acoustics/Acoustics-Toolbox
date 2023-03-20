function plotsr( Pos )

% Plot the source and receiver positions
% usage: plotsr( Pos )
% where Pos is the standard structure containing positions
% e.g. plotray( Pos )
%
% MBP Nov. 2008


% set( gca, 'YDir', 'Reverse' )   % plot with depth-axis positive down
% 
% xlabel( 'Range (km)' )
% xlabel( 'Range (m)' )
% ylabel( 'Depth (m)' )
% title( TITLE )
%figure

hold on

Pos.Nrd = length( Pos.r.z );
Pos.Nrr = length( Pos.r.r );

plot( 1000 * Pos.r.r(  1  ) * ones( Pos.Nsd, 1 ), Pos.s.z, 'ko', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'w', ...
        'MarkerSize', 20 )    % source depths
plot( 1000 * Pos.r.r( end ) * ones( Pos.Nrd, 1 ), Pos.r.z, 'k<', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'w', ...
        'MarkerSize', 5 )    % receiver depths
plot( 1000 * Pos.r.r, Pos.r.z( 1 ) * ones( Pos.Nrr, 1 ),   'kv', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', 'w', ...
        'MarkerSize', 5 )    % receiver ranges

hold off
zoom on
