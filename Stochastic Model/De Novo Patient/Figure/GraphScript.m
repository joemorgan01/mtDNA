%% Graph Script
subplot(2,2,3)
plot(t,100*Sim8./1000,'ko-','MarkerFaceColor','k')
hold on
plot(t,100*Sim9./1000,'ks-','MarkerFaceColor','k')
hold on
plot(t,100*Sim10./1000,'kd-','MarkerFaceColor','k')
hold on
xlabel('Time (years)','FontSize',16);
ylabel('%^{Sim} COX neg cells','FontSize',16)
lgd=legend('10^{-8}','10^{-9}','10^{-10}','Location','northeast');
title(lgd,'Mutation Rate')
axis([0 81 0 inf]);
grid on
set(gca, 'YScale', 'log')