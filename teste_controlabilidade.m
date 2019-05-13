function [valor] = teste_controlabilidade(A,B,Y1,p)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
    valor=1; %supondo que sistema parcialmente controlável
    for i=1:p %varedura para teste de controlabilidade
        teste=Y1(1:end,i)'*B;
        testereal=real(teste);
        testeimag=imag(teste);
        if testereal==0, %realizando operacao teo.autovetor
            valor=0; %caso resultado nulo, pacial nao controlável
        end %finalizando if sum(Y1(i,1:end)'*B)==0
        if testeimag==0, %realizando operacao teo.autovetor
            valor=0; %caso resultado nulo, pacial nao controlável
        end %finalizando if sum(Y1(i,1:end)'*B)==0
    end %finalizando for i=size(1:p)
end

