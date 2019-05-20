function [] = mapeamento_d_est_secao(eigA,eigAfc,alfa,beta,titulo)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    min_y=min(min(imag([eigA,eigAfc]))); %menor parte imag detectada
    max_y=max(max(imag([eigA,eigAfc]))); %maior parte imag detectada
    imagsec=min_y:((max_y-min_y)/99):max_y; %coordenadas imaginarias
    realsec1=-alfa*ones(size(imagsec)); %coordenas reais com alfa
    realsec2=-beta*ones(size(imagsec)); %coordenas reais com alfa
    plot(real(eigA),imag(eigA),'ko',...
         real(eigAfc),imag(eigAfc),'rx',...
         realsec1,imagsec,'b-',...
         realsec2,imagsec,'b-',...
         'MarkerSize',12,...
         'LineWidth',2)
     set(gca,'FontSize',12)
     text=strcat(titulo,' \alpha=',num2str(alfa),' \beta=',num2str(beta));
     title(text,'fontweight','bold')
     xlabel('Real(s)','fontweight','bold')
     ylabel('Imag(s)','fontweight','bold')
     grid on
     legend('Malha Aberta','Malha Fechada','D-regiao')
end

