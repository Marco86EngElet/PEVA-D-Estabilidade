function [] = mapeamento_d_est_disco(eigA,eigAfc,q,r,titulo)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    angulo=0:2*pi/99:2*pi;
    realdisco=r*cos(angulo)-q;
    imagdisco=r*sin(angulo);
     plot(real(eigA),imag(eigA),'ko',...
         real(eigAfc),imag(eigAfc),'rx',...
         realdisco,imagdisco,'b-',...
         'MarkerSize',12,...
         'LineWidth',2)
     set(gca,'FontSize',12)
     text=strcat(titulo,' q=',num2str(q),' r=',num2str(r));
     title(text,'fontweight','bold')
     xlabel('Real(s)','fontweight','bold')
     ylabel('Imag(s)','fontweight','bold')
     grid on
     legend('Malha Aberta','Malha Fechada','D-regiao')
end

