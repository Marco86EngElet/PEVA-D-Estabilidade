function [Lambda_1_til,Y1_til,Q] = complexo_para_real(Lambda_1,Y1,YpA,p)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
        l=0; %iniciando contagem de modos complexos
        for i=1:p %varendo todos os elementos de YpA
            if imag(YpA)~=0 %Caso elemento de YpA tem parte imaginária
                l=l+1; %Conta elemento imaginário
            end %finalizando if imag(YpA)~=0
        end %finalizando for i=1:p 
        l=l/2; %contando pares complexo
        Qj=[1,1;1i,-1i]/sqrt(2); %Matriz Qj
        if l>0 %se há pares complexos
            Q=[]; %Primeiro elemento diagonal
            for i=1:l %vareduras de modos complexo
                Q=blkdiag(Q,Qj); %Construindo parte complexa de Q
            end %finalizando for i=2:l
            if 2*l<p %se quantidade de modos complexo menor que p 
                Q=blkdiag(Q,eye(p-2*l)); %Construindo parte real de Q
            end %finalizando if 2*l<p 
        else %se não há pares complexos
            Q=eye(p); %Construindo parte real de Q
        end %finalizando if l>0
        Y1_til=Y1*Q'/sqrt(2); %Calculando Ytil1 formula 3.9
        Lambda_1_til=Q*Lambda_1*Q'; %Calculando Lambda_til
end

