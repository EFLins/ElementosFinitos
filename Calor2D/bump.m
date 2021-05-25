function r = bump(a,b,n,factor)
  % Use a bump factor to contruct a  stretched space
  % Use is similar to LINSPACE, but nodes are group according to factor
  % if factor > 1 bump nodes to center
  %    factor < 1 bump nodes to edges
  %    factor < 0 & > -1 bump to one edge
  %    factor < -1 & > -inf bump to other  one edge

  if nargin < 3
    n = 100;
  end
  if nargin < 4
    factor = 1;
  end
  if factor > 0 
    x = abs(linspace(-1,1,n)).^factor;
    i2 = ceil(n/2);
    x = [-x(1:i2) x(i2+1:n)];
    r = x/2 + 0.5;
  else
    x = linspace(0,1,n).^abs(factor);
    r = x;
  end 
  r = a + (b-a)*r;

  return