function [F] = D_EVA_secao(A,B,alfa,beta)
    %Opcoes de otimizacao
    options=zeros(1,5); %Valores padroes de opcoes
    options(2)=60; %Numero maximo de iteracoes permitida 
    options(3)=1000000; %Norma euclidiana de [W P] maxima
    options(4)=1; %Porcentagem para velocidade de termino de codigo 
    options(5)=0; %Mostrar execucao da optimizacao feasp
    Atil=A;
    Btil=B;
    p=size(Atil,1);
    m=size(Btil,2);
        setlmis([]); %Inicializando objeto LMI
        [W,nW,sW]=lmivar(2,[m,p]); % definindo W
        [P,nP,sP]=lmivar(1,[p,1]); % definindo P
 		lmiterm([1,1,1,P],2*alfa,1); % 2*alfa*P 
        lmiterm([1,1,1,P],Atil,1,'s'); % Atil*P+P*Atil'
        lmiterm([1,1,1,W],Btil,1,'s'); % Btil*W+W'*Btil'
        lmiterm([1,2,2,P],-2*beta,1); % -2*beta*P
        lmiterm([1,2,2,P],-Atil,1,'s'); % -Atil*P-P*Atil'
        lmiterm([1,2,2,W],-Btil,1,'s'); % -Btil*W-W'*Btil'
        lmiterm([-2,1,1,P],1,1); %P>0
    LMISYS=getlmis; %Fim da estrutura LMI
    [tmin,xopt]= feasp(LMISYS,options); %Buscando Solucao Factivel para LMI
	s = dec2mat(LMISYS,xopt,P); %Determinando P da Solucao Factivel
	w = dec2mat(LMISYS,xopt,W); %Determinando R da Solucao Factivel
    F=w/s;
end
