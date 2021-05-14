%% Graph Script
subplot(2,2,3)
plot(t,100*Sim10./1000,'ko-','MarkerFaceColor','k')
hold on
plot(t,100*Sim11./1000,'ks-','MarkerFaceColor','k')
hold on
plot(t,100*Sim12./1000,'kd-','MarkerFaceColor','k')
hold on
plot(t,100*Sim13./1000,'kv-','MarkerFaceColor','k')
xlabel('Time (years)','FontSize',16);
ylabel('%^{Sim} COX neg cells','FontSize',16)
lgd=legend('10^{-10}','10^{-11}','10^{-12}','10^{-13}','Location','northeast');
title(lgd,'Mutation Rate')
axis([0 81 0 inf]);
grid on
set(gca, 'YScale', 'log')