function res = F_elem(xy, fgen)
% Calcula o vetor elementar do problema de transferencia de calor
% em um dominio 2D. Funcoes de interpolacao bilineares sao utilizadas
% juntamente com integracao de Gauss Legendre


[x, w] = GL_weights(2);

F = zeros(4,1);
for i = 1:numel(x)
  for j = 1:numel(x)
     Fe = func2D(xy,x(i),x(j),fgen);
     F = F + w(i)*w(j)*Fe;
   end
 end
 
res = F;

return

function result = func2D(xy,r,s,fgen)

psi = 0.25*[   (1-r).*(1-s)
               (1+r).*(1-s)
               (1+r).*(1+s)
               (1-r).*(1+s)];

dpsi =  0.25*[      (-1).*(1-s) (1-r).*(-1)
                    ( 1).*(1-s) (1+r).*(-1) 
                    ( 1).*(1+s) (1+r).*( 1)
                    (-1).*(1+s) (1-r).*( 1)];
jac     = dpsi'*xy;
det_jac = det(jac);

x = psi'*xy;

result = psi.*fgen(x)*det_jac;