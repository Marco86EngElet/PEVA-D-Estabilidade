function [F] = D_EVA_disco(A,B,r,q)
    %Opcoes de otimizacao
    options=zeros(1,5); %Valores padroes de opcoes
    options(2)=1000; %Numero maximo de iteracoes permitida 
    options(3)=10000; %Norma euclidiana de [W P] maxima
    options(4)=50; %Porcentagem para velocidade de termino de codigo 
    options(5)=0; %Mostrar execucao da optimizacao feasp
    Atil=A;
    Btil=B;
    n=size(Atil,1);
    m=size(Btil,2);
        setlmis([]); %Inicializando objeto LMI
        [W,nW,sW]=lmivar(2,[m,n]); % definindo W
        [P,nP,sP]=lmivar(1,[n,1]); % definindo P
 		lmiterm([1,1,1,P],-r,1); % -r*P
        lmiterm([1,1,2,P],q,1); % q*P 
        lmiterm([1,1,2,P],Atil,1,'s'); % Atil*P+P*Atil'
        lmiterm([1,1,2,W],Btil,1,'s'); % Btil*W+W'*Btil'
        %lmiterm([1,2,1,P],q,1); % q*P
        %%lmiterm([1,2,1,P],Atil,1,'s'); % -Atil*P-P*Atil'
        %%lmiterm([1,2,1,W],Btil,1,'s'); % -Btil*W-W'*Btil'
        lmiterm([1,2,2,P],-r,1); % -r*P
        lmiterm([-2,1,1,P],1,1); %P>0
        LMISYS=getlmis; %Fim da estrutura LMI
    [tmin,xopt]= feasp(LMISYS,options); %Buscando Solucao Factivel para LMI
	s = dec2mat(LMISYS,xopt,P); %Determinando P da Solucao Factivel
	w = dec2mat(LMISYS,xopt,W); %Determinando R da Solucao Factivel
    F=w/s;
end
