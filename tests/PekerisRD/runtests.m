% Range-dependent bottom type test
% mbp
global units jkpsflag

units = 'km';


%%

bellhop( 'PekerisRDB' )
figure
plotshd( 'PekerisRDB.shd' )
caxisrev( [ 40 80 ] )


%%
% copyfile( 'field_ri.flp', 'field.flp' )
% kraken(  'DickinsFlatK' )
% plotshd( 'DickinsFlatK.shd', 2, 2, 3 )
% caxisrev( [ 70 120 ] )
% 
% scooter( 'DickinsFlatS' )
% plotshd( 'DickinsFlatS.shd.mat', 2, 2, 4 )
% caxisrev( [ 70 120 ] )

% %%
% % Coupled mode runs for range-dependent case
% 
% % run Dickins case, using the wedge ocean
% 
% % cd DickinsPrecalc
% % runtests
% % cd ..
% % 
% % % do it a second way using a simple marching
% % 
% % cd DickinsMarch
% % runtests
% % cd ..
% % 
% % % Note: can see differences between deltar = 1 and 10 m
% ram
% RAMtoSHD
% plotshd( 'DickinsB.shd', 2, 1, 1 )
% caxisrev( [ 70 120 ] )
% plotbty 'DickinsB'
% 
% plotshd( 'SHDFIL.mat', 2, 1, 2 )
% caxisrev( [ 70 120 ] )
% plotbty 'DickinsB'
