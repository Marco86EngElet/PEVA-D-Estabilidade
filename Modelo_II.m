clc,clear all,close all
%----------------------------------------------------------
%Etapa 1
M=[-7.6593,-5.0754,0.9311;6.2936,3.1457,1.2384;...
    -3.5029,-2.4862,-2.0836]; %matriz de massa modelo 1
D=[-2.0374,9.0183,6.6374;0.3073,4.4470,-7.3132;...
     3.1506 ,-1.9984, -8.7907]; %matriz de amortecimento modelo 1
K=[-8.3151,-3.9655,-8.0925;-6.7220,-9.7664,-7.0697;...
    -3.5156,0.7981,2.6228]; %matriz de rigidez Modelo 1
N=eye(3); %matriz de atuadores Bi modelo 1
%-------------------------------------------------
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
%------------------------------------------------
%Etapa 4
disp('Modos de malha aberta escolhidos para alocacao faixa')
YpA_faixa=[7.9807,0.4408,0.7076], %escolhendo os modos
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Verificando quantidade de modos a serem alocados')
p=size(YpA_faixa,2), %quantidade de modos
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
%--------------------------------------------------
%Etapa 5
disp('Matriz Lambda_1_para alocacao faixa')
Lambda_1=diag(YpA_faixa), %Matriz dos modos alocaveis
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Matriz dos Autovetores dos modos escolhidos para alocacao faixa')
Y1=[Y(1:end,1),Y(1:end,5:6)], %Selecionando autovetores desejados
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
%--------------------------------------------------
%Etapa 6 e 7
disp('Verificando se sistema parcialmente controlavel')
[controlavel] = teste_controlabilidade(A,B,p,Y1)
%--------------------------------------------------
%Etapa 8
[Lambda_1_til,Y1_til,Q] = complexo_para_real(Lambda_1,Y1,YpA_faixa,p);
disp('Matriz de Transformação de autovetores para alocacao faixa')
Q, %Mostrando matriz
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Matriz de autovetores na forma real para alocacao faixa')
Y1_til,
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Matriz de autovaloress na forma real para alocacao faixa')
Lambda_1_til,
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
%--------------------------------------------------
%Etapa 9
    alfa=1; %alfa secao para modelo I
    beta=2.5; %beta secao para modelo I
    %LMI Setor
    teta=90; %angulo setor
    angulo=deg2rad(teta); %convertendo angulos em radianos
    %LMI Disco
    r=2;
    q=2;
%--------------------------------------------------
%Etapa 10
disp('Matriz de retroalimentacao para faixa')
Ftil=D_PEVA_secao(Lambda_1_til,B,Y1_til,alfa,beta) %Alocacao Faixa
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Matriz de retroalimentacao para setor')
Ftils=D_PEVA_setor(Lambda_1_til,B,Y1_til,angulo) %Alocacao Setor
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Matriz de retroalimentacao para ')
Ftild=D_PEVA_disco(Lambda_1_til,B,Y1_til,r,q) %Alocacao Disco
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
    titulo='Modelo II - Faixa';
    mapeamento_d_est_secao(eigA,eigAfc,alfa,beta,titulo)
    subplot(132)
    titulo='Modelo II - Setor';
    mapeamento_d_est_setor(eigA,eigAfcs,teta,titulo)
    subplot(133)
    titulo='Modelo II - Disco';
    mapeamento_d_est_disco(eigA,eigAfcd,q,r,titulo)
%--------------------------------------------------