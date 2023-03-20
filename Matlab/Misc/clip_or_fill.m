
function y = clip_or_fill( x, n )

% clips or zero-fills a row or column vector to a prescribed length
% mbp 5/00
% usage: y = clip_or_fill( x );

nx = length( x );

if nx == n		% nothing to do; quick return
  y = x;
  return
end

% row vector

if size( x, 1 ) == 1
  if nx < n
    y = [ x zeros( 1, n - nx ) ];
  end
  
  if nx > n
    y = x( 1:n );
  end
end

% column vector

if size( x, 2 ) == 1
  if nx < n
    y = [ x ; zeros( n - nx, 1 ) ];
  end
  
  if nx > n
    y = x( 1:n );
  end
end
