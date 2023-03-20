function fixfig(afs,lfs,tfs)

if nargin == 0
    afs = 14;	% Axis font size
    lfs = 16;	% Labels font size
    tfs = 18;	% Title font size
end
if nargin == 1
    lfs = afs;
    tfs = afs;
end

fontweights = 'nbb';	% weights for axis, labels, and title ( n = normal, b = bold)


set(gca,'FontName','Arial')
a=findobj(gcf,'Type','axes');
for k=1:length(a)
   set(a(k),'FontSize',afs,'FontWeight',fontweights(1));
   set(get(a(k),'Xlabel'),'FontSize',lfs,'FontWeight',fontweights(2));
   set(get(a(k),'Ylabel'),'FontSize',lfs,'FontWeight',fontweights(2));
   set(get(a(k),'Zlabel'),'FontSize',lfs,'FontWeight',fontweights(2));
   b=get(a(k),'Title');
   set(b,'FontSize',tfs,'FontWeight',fontweights(3));
   c=findobj(get(a(k),'children'),'type','text');
   set(c,'FontSize',lfs,'FontWeight',fontweights(2));
end

a=findobj(gcf,'Tag','suptitle');
if ~isempty(a)
   b=get(a,'children');
   set(b,'FontSize',tfs,'FontWeight',fontweights(3));
end
return


