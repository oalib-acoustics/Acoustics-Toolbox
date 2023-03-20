% Tabulated Reflection Coefficient test cases
% Note that KRAKEN (vs. KRAKENC) cannot use a tabulated reflection coefficient
% mbp

bounce( 'neggradB' )   % make the tabulated reflection coefficient

%%

% KRAKENC runs
krakenc( 'neggradC_geo' )   % reference solution (full geoacoustic bottom)
plotshd( 'neggradC_geo.shd.mat', 3, 1, 1 )
caxisrev( [ 60 130 ] )

copyfile( 'neggradB.brc', 'neggradC_brc.brc' )  % for KRAKENC input
krakenc( 'neggradC_brc' )   % using the BRC file
plotshd( 'neggradC_brc.shd.mat', 3, 1, 2 )
caxisrev( [ 60 130 ] )

copyfile( 'neggradB.irc', 'neggradC_irc.irc' )  % for KRAKENC input
krakenc( 'neggradC_irc' )   % using the IRC file
plotshd( 'neggradC_irc.shd.mat', 3, 1, 3 )
caxisrev( [ 60 130 ] )
%%
% SCOOTER runs

scooter( 'neggradS_geo' )   % reference solution (full geoacoustic bottom)
plotshd( 'neggradS_geo.shd.mat', 3, 1, 1 )
caxisrev( [ 60 130 ] )

movefile( 'neggradB.brc', 'neggradS_brc.brc' )  % for KRAKENC input
scooter( 'neggradS_brc' )   % using the BRC file
plotshd( 'neggradS_brc.shd.mat', 3, 1, 2 )
caxisrev( [ 60 130 ] )

copyfile( 'neggradB.irc', 'neggradS_irc.irc' )  % for KRAKENC input
scooter( 'neggradS_irc' )   % using the IRC file
plotshd( 'neggradS_irc.shd.mat', 3, 1, 3 )
caxisrev( [ 60 130 ] )
