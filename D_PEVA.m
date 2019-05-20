function [F] = D_PEVA(Lambda1,B,Y1,L,V)
    %Opcoes de otimizacao
    options=zeros(1,5); %Valores padroes de opcoes
    options(2)=1000; %Numero maximo de iteracoes permitida 
    options(3)=10000; %Norma euclidiana de [W P] maxima
    options(4)=50; %Porcentagem para velocidade de termino de codigo 
    options(5)=0; %Mostrar execucao da optimizacao feasp
    Atil=Lambda1;
    Btil=Y1'*B;
    p=size(Atil,1);
    m=size(Btil,2);
    setlmis([]); %Inicializando objeto LMI
        [W,nW,sW]=lmivar(2,[m,p]); % definindo W
        [P,nP,sP]=lmivar(1,[p,1]); % definindo P
 		for i=1:max(size(L)) %Laco Selecionando Linha i
			for j=1:max(size(L)) %Laco Selecionando Coluna j
				if L(i,j)~=0 %Verificando se L(i,j) é nao nulo
					lmiterm([1,i,j,P],L(i,j),1); %Lij*P
                end %Encerrando verficacao se L(i,j) é nao nulo
				if V(i,j)~=0 %Verificando se V(i,j) é nao nulo
                    lmiterm([1,i,j,P],V(i,j)*Atil,1,'s')
                    lmiterm([1,i,j,W],V(i,j)*Btil,1,'s'); %Vij*W
                end %Encerrando verficacao se V(i,j) é nao nulo
            end %Fim Laco Selecionando Coluna j
        end %Fim Laco Selecionando Linha i
        lmiterm([-2,1,1,P],1,1); %P>0
    LMISYS=getlmis; %Fim da estrutura LMI
    [tmin,xopt]= feasp(LMISYS,options); %Buscando Solucao Factivel para LMI
	s = dec2mat(LMISYS,xopt,P); %Determinando P da Solucao Factivel
	w = dec2mat(LMISYS,xopt,W); %Determinando R da Solucao Factivel
    F=w/s;
    F=F*Y1';
end
