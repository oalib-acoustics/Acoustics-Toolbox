% Munk profile test cases
% demonstrating role of increasing number of modes
% Fig. 5.12 in JKPS
% mbp
global units jkpsflag
units    = 'km';

% figure
% krakenc( 'MunkK1501' )
% plotshd( 'MunkK1501.shd', 3, 1, 1 )
% caxisrev( [ 60 100 ] )
% title( 'F = 50 Hz, SD = 100 m' )
% 
% krakenc( 'MunkK1510' )
% plotshd( 'MunkK1510.shd', 3, 1, 2 )
% caxisrev( [ 60 100 ] )
% title( '' )
% 
% krakenc( 'MunkKwb' )
% plotshd( 'MunkKwb.shd', 3, 1, 3 )
% caxisrev( [ 60 100 ] )
% title( '' )

%%

krakenc( 'MunkKwb' )
plotshd( 'MunkKwb.shd.mat', 3, 1, 1 )
caxisrev( [ 60 100 ] )
%title( 'F = 50 Hz, SD = 100 m' )

krakenc( 'MunkKbb' )
plotshd( 'MunkKbb.shd.mat', 3, 1, 2 )
caxisrev( [ 60 100 ] )
%title( '' )

krakenc( 'MunkKleaky' )
plotshd( 'MunkKleaky.shd.mat', 3, 1, 3 )
caxisrev( [ 60 100 ] )
%title( '' )

if ( jkpsflag )
   print -depsc fig5_12.eps
end

scooter( 'MunkS' )
figure
plotshd( 'MunkS.shd.mat' )
caxisrev( [ 60 100 ] )
