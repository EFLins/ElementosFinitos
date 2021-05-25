%% inicializacao
clc;
close all;
clearvars; 

tic

k = 0.2; % condutividade
fgen = @(x) 0; % função geracao
fbound = @(x) 0; % funcao fluxo
 
Lx = 1; % dimensoes do dominio
Ly = 1;

Nx = 101; % numero de nos em cada direcao
Ny = 101;

xe = bump(0,Lx,Nx,0.5); % coordenadas nodais
ye = bump(0,Ly,Ny,-0.5);

%% Cria a malha 2D
% matrix com as coordenadas x,y
[xx, yy] = meshgrid(xe,ye); 

% para ver a malha descomente a linha abaixo
% mesh(xx,yy,ones(size(xx))); view(2); return

% numero de nos total
nNos = numel(xx);   
       
% cria uma lista contendo as coordenadas dos nÃ³s em duas colunas [x y]
xnos = [reshape(xx',nNos,1) reshape(yy',nNos,1)]; 
%xnos = xnos+0.02*(rand(size(xnos))-0.5); 

% cria a matriz de malha, onde cada entrada representa um noh
mm = reshape(1:nNos,numel(xe),numel(ye))'; 

% aqui ficara a estrutura da malha, cada linha contem os nos do elemento em sentido horÃ¡rio 
mmesh = [];
% para cada elemento, pega os 4 nos que o define e agrupa na malha
for j = 1:Nx-1
  for i = 1:Ny-1
     mmesh = [mmesh; mm(i,j) mm(i,j+1) mm(i+1,j+1) mm(i+1,j)];
  end
end

%% Definicao das condicoes de contorno
%numero de elementos
nEle = size(mmesh,1);

% inicializas vetores e matriz
U = zeros(nNos,1);
F = zeros(nNos,1);
K = zeros(nNos,nNos);

%obtem os nos no contorno
u_def = unique([mm(:,1); mm(:,end); mm(end,:)'; mm(1,:)']);

du_def = [];

% define os valores para esses nos
% nos com fluxo zero nao precisam ser definido 
U(mm(:,1))   = 0; % esquerda
U(mm(:,end)) = 0; % direita
U(mm(1,:))   = 0; % base
U(mm(end,:)) = 1; % topo

%% Montagem das matrizes 
% matrizes principais
for ele = 1:nEle
  eGl = mmesh(ele,:);
  xy  = xnos(eGl,:);
  klocal = K_elem(xy);
  flocal = F_elem(xy,fgen);
  K(eGl,eGl) = K(eGl,eGl) + klocal;
  F(eGl) = F(eGl) + flocal;
end

% contribuicao do contorno
for ic = 1:numel(du_def)-1
  eGl = du_def([ic ic+1]);
  xy  = xnos(eGl,:);
  flocal = bound_elem(xy,fbound);
  F(eGl) = F(eGl) + flocal;
end

%% solucao
nos_free = setdiff(1:nNos,u_def);
U(nos_free) = K(nos_free,nos_free)\...
   (F(nos_free)-K(nos_free,u_def)*U(u_def));

toc

%% pos-processamento
UU = reshape(U,Nx,Ny)';

figure;
mesh(xx,yy,UU,'EdgeColor','b','FaceAlpha',0.1);  hold on; view(2);
% quad_display([],mmesh,xnos',U);
contourf(xx,yy,UU,[0:0.1:1],'linestyle','-','Showtext','on'); axis equal
colormap(parula); colorbar;
set(gca(), "fontsize",16);
xlabel("X");
ylabel("Y");
axis square;


figure
sa = load('Ex1_analitic.mat');
plot(yy(xx==Lx/2),UU(xx==Lx/2),'bo-',...
    sa.y(sa.x==Lx/2),sa.T(sa.x==Lx/2),'r-')
xlabel("X")
ylabel("Temperatura");
set(gca(), "fontsize",16)
