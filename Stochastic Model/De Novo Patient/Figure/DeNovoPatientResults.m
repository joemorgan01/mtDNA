close all
load('1kResults');
mutation=100*(Results(:,2,:))./(Results(:,2,:)+Results(:,1,:));

% Need to change the mutation array to be 1000 x 1 x 4160 instead of 
% 4160 x 1 x 1000
Dimensions=size(mutation);

Mutation4KSDensity=NaN(Dimensions(3), Dimensions(2), Dimensions(1));
%                        1000       x        1     x      4160
Densities=NaN(1,100,Dimensions(1));
MutationLoads=NaN(1,100,Dimensions(1));

for i=1:Dimensions(3)
    Mutation4KSDensity(i,:,:)=mutation(:,:,i);
end

for i=1:Dimensions(1)
[Density, MutationLoad]= ksdensity(Mutation4KSDensity(:,:,i),'BandWidth',0.2);
Densities(:,:,i)=Density;
MutationLoads(:,:,i)=MutationLoad;
end

Dens=NaN(4160,100);

for i=1:Dimensions(1)
    Dens(i,:)=Densities(:,:,i);
end
   

 x=linspace(-3, 103, 100);
 y=1:4160;
 Z=Dens;

subplot(2,2,1)
surface(x,y./52,Dens,'FaceAlpha',0.8, 'EdgeColor','none');
colormap jet
grid on
axis([0 100 0.05 80 0 0.7]);
xlabel('Mutation Load (%)','FontSize',14);
ylabel('Time (years)','FontSize',12);
zlabel('Density','FontSize',14);
set(gcf,'color','w')
az = 13;
el = 21;
view([az,el])

clear

load('Sim100');
load('Sim1000');
load('Sim2000');
load('Sim5000');
load('Sim10000');
load('time');
run('RunNoGraphScript');

clear

load('Sim10');
load('Sim9');
load('Sim8');
load('time');
run('GraphScript');

clear

load('y');
run('TimeofFounder');

set(gcf,'position',[681,320,898,657])
print(gcf,'DeNovoPatient.png','-dpng','-r400')
