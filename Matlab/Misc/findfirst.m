function Iloc = findfirst( x )
% findfirst( x)
% finds the position of the first one in the logical vector x
% assumes (!!!) x is monotonic, i.e. x is a block of zeros followed by
% a block on ones
% uses a bisection approach

ileft = 1;
iright = length( x );

while ( ileft + 1 < iright )
  imid = round( ( ileft + iright ) / 2 );
  if x( imid ) == 0
    ileft = imid;
  else
    iright = imid;
  end
end

Iloc = iright;



