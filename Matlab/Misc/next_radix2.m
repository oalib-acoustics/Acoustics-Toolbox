function M = next_radix2( L )
n=1;
N=2;
while L > N
  n = n+1;
  N = 2^n;
end
M = N;
