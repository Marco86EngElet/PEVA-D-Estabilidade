clc,clear all,close all
%----------------------------------------------------------
%Etapa 1
M=[17.6,1.28,2.89; 1.28,0.824,0.413;2.89,0.413,0.725]; %matriz de massa modelo 1
D=[7.66,2.45,2.1;0.23,1.04,0.223;0.60,0.756,0.658]; %matriz de amortecimento modelo 1
K=[121,18.9,15.9;0,2.7,0.145;11.9,3.64,15.5];
N=[1;1;1]; %matriz de atuadores Bi modelo 1
%----------------------------------------------
%Etapa 2
n=size(M,1); %determinando grau de liberdade
m=size(N,2); %determinando ordem de B
A=[zeros(n),eye(n);-inv(M)*K,-inv(M)*D]; %Matrix de estados
B=[zeros(n,m);inv(M)*N]; %matrix B de entradas
%--------------------------------------------------
%Etapa 3
[Vd1,Lambda,Y]=eig(A); %Determinando autovetores e autovalores de A
disp('Autovalores em malha aberta')
diag(Lambda), %mostrrando matriz Lambda
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
%--------------------------------------------------
%Etapa 4
disp('Modos de malha aberta escolhidos para alocacao faixa')
YpA_faixa=[0.0947 + 2.5229i,0.0947 - 2.5229i], %escolhendo os modos
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Verificando quantidade de modos a serem alocados')
p_faixa=size(YpA_faixa,2), %quantidade de modos
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Modos de malha aberta escolhidos para alocacao setor e disco')
YpA_setor=[-0.8848 + 8.4415i,-0.8848 - 8.4415i,YpA_faixa], %escolhendo os modos
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Verificando quantidade de modos a serem alocados')
p_setor=size(YpA_setor,2), %quantidade de modos
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
%--------------------------------------------------
%Etapa 5
disp('Matriz Lambda_1_para alocacao faixa')
Lambda_1_faixa=diag(YpA_faixa), %Matriz dos modos alocaveis
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Matriz dos Autovetores dos modos escolhidos para alocacao faixa')
Y1_faixa=Y(1:end,3:4), %Selecionando autovetores desejados
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Matriz Lambda_1_para alocacao setor')
Lambda_1_setor=diag(YpA_setor), %Matriz dos modos alocaveis
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Matriz dos Autovetores dos modos escolhidos para alocacao setor')
Y1_setor=Y(1:end,1:4), %Selecionando autovetores desejados
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
%--------------------------------------------------
%Etapa 6 e 7
disp('Verificando se sistema parcialmente controlavel faixa')
[controlavel] = teste_controlabilidade(A,B,p_faixa,Y1_faixa)
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Verificando se sistema parcialmente controlavel setor')
[controlavel] = teste_controlabilidade(A,B,p_setor,Y1_setor)
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
%--------------------------------------------------
%Etapa 8
[Lambda_1_til_faixa,Y1_til_faixa,Q_faixa] = ...
    complexo_para_real(Lambda_1_faixa,Y1_faixa,YpA_faixa,p_faixa);
disp('Matriz de Transformação de autovetores para alocacao faixa')
Q_faixa, %Mostrando matriz
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Matriz de autovetores na forma real para alocacao faixa')
Y1_til_faixa,
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Matriz de autovaloress na forma real para alocacao faixa')
Lambda_1_til_faixa,
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
[Lambda_1_til_setor,Y1_til_setor,Q_setor] = ...
    complexo_para_real(Lambda_1_setor,Y1_setor,YpA_setor,p_setor);
disp('Matriz de Transformação de autovetores para alocacao setor')
Q_setor, %Mostrando matriz
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Matriz de autovetores na forma real para alocacao setor')
Y1_til_setor,
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Matriz de autovaloress na forma real para alocacao setor')
Lambda_1_til_setor,
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
%--------------------------------------------------
%Etapa 9
    alfa=0.4; %alfa secao para modelo I
    beta=0.6; %beta secao para modelo I
    %LMI Setor
    teta=110; %angulo setor
    angulo=deg2rad(teta); %convertendo angulos em radianos
    %LMI Disco
    r=0.5; %raio
    q=1.5; %centro
%--------------------------------------------------
%Etapa 10
disp('Matriz de retroalimentacao para faixa')
Ftil=D_PEVA_secao(Lambda_1_til_faixa,B,Y1_til_faixa,alfa,beta) %Alocacao Faixa
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Matriz de retroalimentacao para setor')
Ftils=D_PEVA_setor(Lambda_1_til_setor,B,Y1_til_setor,angulo) %Alocacao Setor
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Matriz de retroalimentacao para disco ')
Ftild=D_PEVA_disco(Lambda_1_til_setor,B,Y1_til_setor,r,q) %Alocacao Disco
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
%--------------------------------------------------
%Etapa 11
eigA=cplxpair(polyeig(K,D,M));
eigAfc=cplxpair(polyeig(K-N*Ftil(1:m,1:n),...
    D-N*Ftil(1:m,(n+1):2*n),M));
eigAfcd=cplxpair(polyeig(K-N*Ftild(1:m,1:n),...
                      D-N*Ftild(1:m,(n+1):2*n),M));
eigAfcs=cplxpair(polyeig(K-N*Ftils(1:m,1:n),...
                      D-N*Ftils(1:m,(n+1):2*n),M));
%--------------------------------------------------
%Etapa 12
fprintf('  Modos Malha Aberta\t\t\t\t\tModos Malha Fechada')
fprintf('\n\t\t\t\t\t\t\tSecao\t\t\t  Setor\t\t\t\tDisco\n')
disp([eigA,eigAfc,eigAfcs,eigAfcd])
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
%--------------------------------------------------
%Etapa 13
figure
    set(gcf,'color','w')
    subplot(131)
    titulo='Modelo III - Faixa';
    mapeamento_d_est_secao(eigA,eigAfc,alfa,beta,titulo)
    subplot(132)
    titulo='Modelo III - Setor';
    mapeamento_d_est_setor(eigA,eigAfcs,teta,titulo)
    subplot(133)
    titulo='Modelo III - Disco';
    mapeamento_d_est_disco(eigA,eigAfcd,q,r,titulo)
%--------------------------------------------------
F1=Ftild;
F2=D_EVA_disco(A,B,r,q);
k=1.025