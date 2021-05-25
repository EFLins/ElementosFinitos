function res = K_elem(xy)

% carrega os pontos e abcissas da integracao de Gauss
[x, w] = GL_weights(2);

K = zeros(4);

for i = 1:numel(x)
  for j = 1:numel(x)
     Ke = func2D(xy,x(i),x(j));
     K = K + w(i)*w(j)*Ke;
   end
 end
 
res = K;

function result = func2D(xy,r,s)

psi = 0.25*[(1-r).*(1-s)
            (1+r).*(1-s)
            (1+r).*(1+s)
            (1-r).*(1+s)];

dpsi = 0.25*[(-1).*(1-s) (1-r).*(-1)
             ( 1).*(1-s) (1+r).*(-1) 
             ( 1).*(1+s) (1+r).*( 1)
             (-1).*(1+s) (1-r).*( 1)];
jac     = dpsi'*xy;
det_jac = det(jac);
%ijac    = inv(jac);
ijac    = [jac(2,2) -jac(1,2); -jac(2,1) jac(1,1)]/det_jac;
B       = ijac*dpsi';
result = B'*det_jac*B;