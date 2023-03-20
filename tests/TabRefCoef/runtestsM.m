% Tabulated Reflection Coefficient test cases
% mbp

bounce( 'neggradB' )   % make the tabulated reflection coefficient

%krakenc( 'neggradC_geo' )   % reference solution (full geoacoustic bottom)
%plotshd( 'neggradC_geo.shd', 3, 1, 1 )
%caxisrev( [ 60 130 ] )

%copyfile( 'neggradB.brc', 'neggradC_brc.brc' )  % for KRAKENC input
%krakenc( 'neggradC_brc' )   % using the BRC file
%plotshd( 'neggradC_brc.shd', 3, 1, 2 )
%caxisrev( [ 60 130 ] )

%copyfile( 'neggradB.irc', 'neggradC_irc.irc' )  % for KRAKENC input
%krakenc( 'neggradC_irc' )   % using the IRC file
%plotshd( 'neggradC_irc.shd', 3, 1, 3 )
%caxisrev( [ 60 130 ] )

scooterM( 'neggradS_geo' )   % reference solution (full geoacoustic bottom)
plotshd( 'neggradS_geo.shd.mat', 2, 1, 1 )
caxisrev( [ 60 130 ] )

movefile( 'neggradB.brc', 'neggradS_brc.brc' )  % for KRAKENC input
scooterM( 'neggradS_brc' )   % using the BRC file
plotshd( 'neggradS_brc.shd.mat', 2, 1, 2 )
caxisrev( [ 60 130 ] )

%copyfile( 'neggradB.irc', 'neggradS_irc.irc' )  % for KRAKENC input
%scooterM( 'neggradS_irc' )   % using the IRC file
%plotshd( 'neggradS_irc.shd.mat', 3, 1, 3 )
%caxisrev( [ 60 130 ] )