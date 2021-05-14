%% Graph Script
subplot(2,2,2)
plot(t,100*Sim100./100,'ko-','MarkerFaceColor','k')
hold on
plot(t,100*Sim1000./1000,'ks-','MarkerFaceColor','k')
hold on
plot(t,100*Sim2000./2000,'kd-','MarkerFaceColor','k')
hold on
plot(t,100*Sim5000./5000,'kv-','MarkerFaceColor','k')
hold on
plot(t,100*Sim10000./10000,'kp-','MarkerFaceColor','k')
xlabel('Time (years)','FontSize',16);
ylabel('%^{Sim} COX neg cells','FontSize',16)
legend('N=100','N=1000','N=2000','N=5000','N=10000','Location','northwest');
axis([0 81 0 inf]);
grid on