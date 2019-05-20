figure
    set(gcf,'color','w')
    subplot(131)
    titulo='h_1(t)';
    plot(ht_PEVA(1:end,1),ht_PEVA(1:end,2),'b-',...
         ht_EVA(1:end,1),ht_EVA(1:end,2),'r-','LineWidth',2)
     set(gca,'FontSize',12)
     text=titulo; title(text,'fontweight','bold'),
     xlabel('t','fontweight','bold'), ylabel('h_1','fontweight','bold')
     grid on, legend('PEVA','EVA')
    subplot(132)
    titulo='h_2(t)';
    plot(ht_PEVA(1:end,1),ht_PEVA(1:end,3),'b-',...
         ht_EVA(1:end,1),ht_EVA(1:end,3),'r-','LineWidth',2)
     set(gca,'FontSize',12), text=titulo;
     title(text,'fontweight','bold'), xlabel('t','fontweight','bold')
     ylabel('h_2','fontweight','bold'), grid on
     legend('PEVA','EVA')
    subplot(133)
    titulo='h_3(t)';
    plot(ht_PEVA(1:end,1),ht_PEVA(1:end,4),'b-',...
         ht_EVA(1:end,1),ht_EVA(1:end,4),'r-','LineWidth',2)
     set(gca,'FontSize',12)
     text=titulo;
     title(text,'fontweight','bold')
     xlabel('t','fontweight','bold')
     ylabel('h_3','fontweight','bold')
     grid on
     legend('PEVA','EVA')