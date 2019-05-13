function [] = mapeamento_d_est_setor(eigA,eigAfc,teta,titulo)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    min_x=min(min(real([eigA,eigAfc]))); %menor parte real detectada
    min_y=min(min(imag([eigA,eigAfc]))); %menor parte imag detectada
    max_x=max(max(real([eigA,eigAfc]))); %maior parte real detectada
    max_y=max(max(imag([eigA,eigAfc]))); %maior parte imag detectada
    realset=0:min_x/99:min_x; %coordenadas x para o setor
    imagset1=realset*tand(teta/2);
    imagset2=-realset*tand(teta/2);
    plot(real(eigA),imag(eigA),'ko',...
         real(eigAfc),imag(eigAfc),'rx',...
         realset,imagset1,'b-',...
         realset,imagset2,'b-',...
         'MarkerSize',12,...
         'LineWidth',2)
     set(gca,'FontSize',12)
     text=strcat(titulo,' \theta=',num2str(teta),'\circ');
     title(text,'fontweight','bold')
     xlabel('Real(s)','fontweight','bold')
     ylabel('Imag(s)','fontweight','bold')
     grid on
    legend('Malha Aberta','Malha Fechada','D-regiao')

end

