function x = tridiag2( a, b, c, y )
%      x = tridiag2( a, b, c, y )
%
%  Solve the  N x N  tridiagonal system for x:
%
%  [ b(1)  c(1)                                  ] [  x(1)  ]   [  y(1)  ] 
%  [ a(1)  b(2)  c(2)                            ] [  x(2)  ]   [  y(2)  ] 
%  [       a(2)  b(3)  c(3)                      ] [        ]   [        ] 
%  [            ...   ...   ...                  ] [  ...   ] = [  ...   ]
%  [                    ...    ...    ...        ] [        ]   [        ]
%  [                        a(N-2) b(N-1) c(N-1) ] [ x(N-1) ]   [ y(N-1) ]
%  [                               a(N-1)  b(N)  ] [  x(N)  ]   [  y(N)  ]
%
%  y must be a vector (row or column); N is determined from its length.
%  a, b, c must be vectors of lengths at least N-1, N, and N-1 respectively,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Check that the input arrays have acceptable sizes

if min(size(y))~=1 | min(size(a))~=1 | min(size(b))~=1 | min(size(c))~=1
   error('a, b, c, y must be vectors');
end
N = length(y);
if length(a)<N-1 | length(b)<N | length(c)<N-1
   error('a, b, c must be vectors of length at least N-1, N, N-1');
end

isrow = ( size(y,2)>size(y,1) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Solve the problem using Matlab's sparse matrix routines

% Make them all be column vectors if they are rows
a=a(:);  b=b(:);  c=c(:);  y=y(:);

%  Form the tridiagonal sparse matrix
A = spdiags( [ [ a(1:N-1); 0 ] b(1:N) [ 0; c(1:N-1) ] ], -1:1, N, N );

%  Solve the system using Matlab's sparse solver
x = A \ y;

%  Make the returned array x have the same orientation as input y
if isrow;   x = x';  end;
