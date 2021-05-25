function res = bound_elem(xy, fbound)
% Calcula o vetor fluxo elementar do problema de transferencia de calor
% em um dominio 2D. Funcoes de interpolacao lineares sao utilizadas
% juntamente com integracao de Gauss Legendre

[x, w] = GL_weights(2);

F = zeros(2,1);

for i = 1:numel(x)
  Fe = func1D(xy,x(i),fbound);
  F = F + w(i)*Fe;
end
 
res = F;

return

function result = func1D(xy2D,r,fbound)

% o tamanho do elemento eh a diferenca 
% entre as coordenadas dos nos
h = norm(diff(xy2D));

psi = 0.5*[(1-r)
           (1+r)];
                    
det_jac = h/2;
x = psi'*xy2D; % obtem os valores (x,y) no espaco fisico
result = psi*fbound(x)*det_jac;