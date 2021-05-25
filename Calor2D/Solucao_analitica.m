function Solucao_analitica(ishow)
% Solucao do problema de transferencia de calor 2D
% utiliza separacao de variaveis. Os resultados sao salvos em um
% arquivo .mat que pode ser lido por outra funcao.
% a variavel ishow deve ser 1 para ocultar 
% os graficos na tela

if nargin == 0
 ishow = 0;
 end
 

L = 1;
W = 1;

[x, y] = meshgrid(...
    linspace(0,L,101),...
    linspace(0,W,101));

T = zeros(size(x));

for n = 1:200
   a = ((-1)^(n+1)+1)/n;
   b = sin(n*pi*x/L);
   c = sinh(n*pi*y/L)/sinh(n*pi*W/L);
   
   T = T + 2/pi*a*b.*c;
   
end

if ishow ~= 1
	figure;
	contourf(x,y,T,0:0.1:1,'ShowText','on');
	colormap(jet); colorbar;
	set(gca(),"fontsize",16)
	xlabel("X")
	ylabel("Y")
    title("Temperatura")
	axis square

	figure;
	plot(y(x==L/2),T(x==L/2),'bo-')
	xlabel("X")
	ylabel("Temperatura")
	set(gca(),"fontsize",16)
end
save('Ex1_analitic.mat','x','y','T');

return
