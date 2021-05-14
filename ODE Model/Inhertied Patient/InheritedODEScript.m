clear all
close all

year=365*24*3600;

StopTime=80*year;
MoleculeNo=200;
MutantPercent=0.65;

InitialConditions=[MoleculeNo*(1-MutantPercent) MoleculeNo*MutantPercent];
[t,C]=ode45(@InheritedODEModel,[0 StopTime],InitialConditions);

figure(1)
plot(t./year,C);
axis([0 StopTime/year 0 MoleculeNo]);
xlabel('Time / years');
ylabel('Copy Number');
legend('Wild-Type', 'Mutant','Location', 'northeast');
grid on
saveas(figure(1),'Inherited ODE - Component Copy Number','jpeg');

figure(2)
mutation=100*C(:,2)./(C(:,2)+C(:,1));
plot(t./year,mutation);
axis([0 StopTime/year 0 100]);
xlabel('Time / years');
ylabel('Mutation %');
title('Mutation %')
grid on
saveas(figure(2),'Inherited ODE - Mutation %','jpeg');

figure(3)
totalcopyno=C(:,2)+C(:,1);
plot(t./year, totalcopyno)
axis([0 StopTime/year 0 1.1*MoleculeNo]);
xlabel('Time / years');
ylabel('Copy Number');
title('Total Copy Number');
grid on
saveas(figure(3),'Inherited ODE - Total Copy Number','jpeg');