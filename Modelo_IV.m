clc,clear all,close all
%----------------------------------------------------------
%Etapa 1
M=diag([0.7177,1.2660,1.3918],1);
M=M+diag([0.4757,0.9676],2);
M(1,4)=0.4311;
M=M+M';
M=M+diag([1.4685,2.6938,2.7061,2.1876]); %matriz de massa
K=diag([0.0076,-0.0101,-0.2564],1);
K=K+diag([-0.1359,-0.0493],2);
K(1,4)=-0.7290;
K=K+K';
K=K+diag([1.7824,1.0287,2.8360,1.9130]); %matriz de rigidez
D=diag([1.2695,0.9144,0.8310],1);
D=D+diag([0.7967,0.7325],2);
D(1,4)=0.8160;
D=D+D';
D=D+diag([1.3525,1.3274,0.9456,1.1536]); %matriz de amortecimento
N=[1;2;1;-1];
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
YpA_faixa=[-0.0861 + 1.6242i,-0.0861 - 1.6242i], %escolhendo os modos
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Verificando quantidade de modos a serem alocados')
p_faixa=size(YpA_faixa,2), %quantidade de modos
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Modos de malha aberta escolhidos para alocacao setor e disco')
YpA_setor=[YpA_faixa,-0.1748+1.1922i,-0.1748-1.1922i,...
            -0.1022 + 0.8876i,-0.1022-0.8876i], %escolhendo os modos
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
Y1_faixa=Y(1:end,1:2), %Selecionando autovetores desejados
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Matriz Lambda_1_para alocacao setor')
Lambda_1_setor=diag(YpA_setor), %Matriz dos modos alocaveis
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Matriz dos Autovetores dos modos escolhidos para alocacao setor')
Y1_setor=[Y(1:end,1:4),Y(1:end,7:8)], %Selecionando autovetores desejados
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
%--------------------------------------------------
%Etapa 6 e 7
disp('Verificando se sistema parcialmente controlavel Faixa')
[controlavel] = teste_controlabilidade(A,B,p_faixa,Y1_faixa)
disp('%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%')
disp('Verificando se sistema parcialmente controlavel Setor')
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
    alfa=0.3; %alfa secao para modelo I
    beta=0.4; %beta secao para modelo I
    %LMI Setor
    teta=120; %angulo setor
    angulo=deg2rad(teta); %convertendo angulos em radianos
    %LMI Disco
    r=1.5; %raio
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
    titulo='Modelo IV - Faixa';
    mapeamento_d_est_secao(eigA,eigAfc,alfa,beta,titulo)
    subplot(132)
    titulo='Modelo IV - Setor';
    mapeamento_d_est_setor(eigA,eigAfcs,teta,titulo)
    subplot(133)
    titulo='Modelo IV - Disco';
    mapeamento_d_est_disco(eigA,eigAfcd,q,r,titulo)