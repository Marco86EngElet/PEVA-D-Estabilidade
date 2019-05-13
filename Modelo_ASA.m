clc,clear all,close all
%----------------------------------------------------------
%Modelo IV Entradas
M1=diag([0.7177,1.2660,1.3918],1);
M1=M1+diag([0.4757,0.9676],2);
M1(1,4)=0.4311;
M1=M1+M1';
M1=M1+diag([1.4685,2.6938,2.7061,2.1876]); %matriz de massa
D1=[7.66,2.45,2.1;0.23,1.04,0.223;0.60,0.756,0.658]; %matriz de amortecimento modelo 1
K1=diag([0.0076,-0.0101,-0.2564],1);
K1=K1+diag([-0.1359,-0.0493],2);
K1(1,4)=-0.7290;
K1=K1+K1';
K1=K1+diag([1.7824,1.0287,2.8360,1.9130]); %matriz de rigidez
D1=diag([1.2695,0.9144,0.8310],1);
D1=D1+diag([0.7967,0.7325],2);
D1(1,4)=0.8160;
D1=D1+D1';
D1=D1+diag([1.3525,1.3274,0.9456,1.1536]); %matriz de amortecimento
B1=[1;2;1;-1];
%------------------------------------------------------
    r=0.5; %raio
    q=1.5; %centro
%------------------------------------------------------
n=size(M1,1); %determinando grau de liberdade
m=size(B1,2); %determinando ordem de B
A=[zeros(n),eye(n);-inv(M1)*K1,-inv(M1)*D1]; %Matrix de estados
B=[zeros(n,m);inv(M1)*B1]; %matrix B de entradas
[Vd1,Lambda,Y]=eig(A); %Determinando autovetores e autovalores de A
%------------------------------------------------------
disp('Autovalores em malha aberta')
diag(Lambda), %mostrrando matriz Lambda
%---------------------------------------------------------
disp('Modos de malha aberta escolhidos para alocacao faixa')
YpA_faixa=[-0.0861 + 1.6242i,-0.0861 - 1.6242i], %escolhendo os modos
disp('Modos de malha aberta escolhidos para alocacao setor e disco')
YpA_setor=[YpA_faixa,-0.1748+1.1922i,-0.1748-1.1922i,...
            -0.1022 + 0.8876i,-0.1022-0.8876i], %escolhendo os modos
disp('Verificando quantidade de modos a serem alocados')
p_setor=size(YpA_setor,2), %quantidade de modos
disp('Matriz Lambda_1_para alocacao setor')
Lambda_1_setor=diag(YpA_setor), %Matriz dos modos alocaveis
disp('Matriz dos Autovetores dos modos escolhidos para alocacao setor')
Y1_setor=[Y(1:end,1:4),Y(1:end,7:8)], %Selecionando autovetores desejados
[Lambda_1_til_setor,Y1_til_setor,Q_setor] = ...
    complexo_para_real(Lambda_1_setor,Y1_setor,YpA_setor,p_setor);
disp('Matriz de Transformação de autovetores para alocacao setor')
Q_setor, %Mostrando matriz
disp('Matriz de autovetores na forma real para alocacao setor')
Y1_til_setor,
disp('Matriz de autovaloress na forma real para alocacao setor')
Lambda_1_til_setor,

%---------------------------------------------------------
disp('Matriz de retroalimentacao para disco via EVA ')
F_EVA=D_EVA_disco(A,B,r,q) %Alocacao Disco
%---------------------------------------------------------
disp('Matriz de retroalimentacao para disco via PEVA ')
F_PEVA=D_PEVA_disco(Lambda_1_til_setor,B,Y1_til_setor,r,q) %Alocacao Disco
