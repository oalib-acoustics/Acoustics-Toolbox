% plots a single TL surface in dB
% compare results to Fig. 6.8, JKPS on p. 376

figure
plotshd( 'wedge.shd' )
caxisrev( [ 40 80 ] )
%set(1,'PaperPosition', [ 0.25 0.25 5.0 3.5 ] )
%print -deps tl.ps

figure
plottlr( 'wedge.shd', 150 )
axis( [ 0 4 40 90 ] );

figure
plottlr( 'wedge.shd',  30 )
axis( [ 0 4 40 90 ] );

%set(2,'PaperPosition', [ 0.25 0.25 5.0 3.5 ] )
%print -deps tlslice.ps
