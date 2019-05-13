clc,clear all,close all
%----------------------------------------------------------
%Etapa 1
M=eye(211);
D=-0.1*diag(ones(1,210),1);
D=D+D';
D=D+diag([0.4,0.5*ones(1,209),0.4]);
K=-100*diag([ones(1,209),0],1);
K=K-100*diag(ones(1,210),-1);
K=K+200*diag([0.5,ones(1,209),0.5]);
N=eye(211,1);
%--------------------------------------------------------
%Etapa 2
n=size(M,1); %determinando grau de liberdade
m=size(N,2); %determinando ordem de B
A=[zeros(n),eye(n);-inv(M)*K,-inv(M)*D]; %Matrix de estados
B=[zeros(n,m);inv(M)*N]; %matrix B de entradas
%--------------------------------------------------
%Etapa 3
[Vd1,Lambda,Y]=eig(A); %Determinando autovetores e autovalores de A
disp('Autovalores em malha aberta (Omitido)')
diag(Lambda); %mostrrando matriz Lambda
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%');
%--------------------------------------------------
%Etapa 4
disp('Modos de malha aberta escolhidos para alocacao faixa')
resp=diag(Lambda);
[row,col]=find(resp>-0.02)
YpA_faixa=Lambda(row,row), %escolhendo os modos
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Verificando quantidade de modos a serem alocados')
p=size(YpA_faixa,2), %quantidade de modos
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
%--------------------------------------------------
%Etapa 5
disp('Matriz Lambda_1_para alocacao faixa')
Lambda_1=diag(YpA_faixa), %Matriz dos modos alocaveis
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Matriz dos Autovetores dos modos escolhidos (Omitido)')
Y1=Y(1:end,row); %Selecionando autovetores desejados
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
%----------------------------------------
%Etapa 6 e 7
disp('Verificando se sistema parcialmente controlavel')
[controlavel] = teste_controlabilidade(A,B,p,Y1)
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
%---------------------------------------- 
%Etapa 8
[Lambda_1_til,Y1_til,Q] = complexo_para_real(Lambda_1,Y1,YpA_faixa,p);
disp('Matriz de Transformação de autovetores para alocacao')
Q, %Mostrando matriz
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Matriz de autovetores na forma real para alocacao (Omitido)')
Y1_til;
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Matriz de autovaloress na forma real para alocacao')
Lambda_1_til,
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
%----------------------------------------
%Etapa 9
    %LMI Faixa
    alfa=0.1; %alfa secao para modelo 
    beta=0.2; %beta secao para modelo 
%---------------------------------------------------------
%Etapa 10
disp('Matriz de retroalimentacao para faixa (Omitido)')
Ftil=D_PEVA_secao(Lambda_1_til,B,Y1_til,alfa,beta); %Alocacao Faixa
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
%------------------------------------------------
%Etapa 11
eigA=cplxpair(polyeig(K,D,M));
eigAfc=cplxpair(polyeig(K-N*Ftil(1:m,1:n),...
    D-N*Ftil(1:m,(n+1):2*n),M));
%------------------------------------------------
%Etapa 12
fprintf('  Modos Malha Aberta\t\tModos Malha Fechada')
fprintf('\n\t\t\t\t\t\t\tSecao\n')
disp([eigA(end-5:end),eigAfc(end-5:end)])
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
%-------------------------------------------------
%Etapa 13
figure
    set(gcf,'color','w')
    titulo='Modelo VI - Faixa';
    mapeamento_d_est_secao(eigA,eigAfc,alfa,beta,titulo)