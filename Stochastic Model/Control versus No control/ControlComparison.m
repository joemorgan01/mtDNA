clear all;
close all;
year=365*24*3600;
day=24*3600;
rxnT=80*year;
%rxnT=30*day;

RunNo=1000;

sampletime=linspace(0,rxnT,52*rxnT/year);

% Initialise Array Sizes
Yplot200=NaN(2, length(sampletime));
tplot=NaN;
Results=NaN(length(sampletime),2 , RunNo);
TResults=NaN(1,length(sampletime), RunNo);
COXDefs=NaN(length(sampletime),1,RunNo);
Inlet=NaN(1,2,RunNo);
Rates=NaN(1,5,RunNo);
Mol=NaN(1,1,RunNo);

% Determine simulation condition by sampling from probability distribution
% funcitons
for i=1:RunNo
    MoleculeNo=InitialCopyNo();
%MoleculeNo=14;
    Y=InitialConditionsv2(MoleculeNo);
    C0=RatesWithUncertaincy();
    Rates(:,:,i)=C0;
    Inlet(:,:,i)=Y;
    Mol(:,:,i)=MoleculeNo;
end
disp('Values Made');
%Run Simulation
for i=1:RunNo
    InitialCondition=Inlet(:,:,i);
    Rate=Rates(:,:,i);
    [tplot,Yplot200,COXDef] = mtDNAModelv2noctrl([InitialCondition(1) InitialCondition(2)], sampletime, Rate);
    Results(:,:,i) = Yplot200;
    TResults(:,:,i)=tplot;
    COXDefs(:,:,i)=COXDef;
    i=i
end

copynumber=(Results(:,2,:)+Results(:,1,:));

% Plot Results
subplot(1,2,1)
avcopyno=mean(copynumber,3);
copyno95=prctile(copynumber,95,3);
copyno5=prctile(copynumber,5,3);
stairs((sampletime./year),avcopyno,'k-','LineWidth',3);
hold on
stairs((sampletime./year),copyno95,'k:','LineWidth',0.5);
hold on
stairs((sampletime./year),copyno5,'k:','LineWidth',0.5);
%legend('Mean','95^{th} Percentile','5^{th} Percentile','Location','northwest','FontSize',14);
grid on
ylabel('Copy Number','FontSize',16);
xlabel('Time (years)','FontSize',16);
title('No Control','FontSize',16);
axis([0 inf 0 inf]);

%Run Simulation
for i=1:RunNo
    InitialCondition=Inlet(:,:,i);
    Rate=Rates(:,:,i);
    [tplot,Yplot200,COXDef] = mtDNAModelv2wControl([InitialCondition(1) InitialCondition(2)], sampletime, Rate);
    Results(:,:,i) = Yplot200;
    TResults(:,:,i)=tplot;
    COXDefs(:,:,i)=COXDef;
    i=i
end

copynumber=(Results(:,2,:)+Results(:,1,:));

subplot(1,2,2)
avcopyno=mean(copynumber,3);
copyno95=prctile(copynumber,95,3);
copyno5=prctile(copynumber,5,3);
stairs((sampletime./year),avcopyno,'k-','LineWidth',3);
hold on
stairs((sampletime./year),copyno95,'k:','LineWidth',0.5);
hold on
stairs((sampletime./year),copyno5,'k:','LineWidth',0.5);
%legend('Mean','95^{th} Percentile','5^{th} Percentile','Location','northwest','FontSize',14);
grid on
ylabel('Copy Number','FontSize',16);
xlabel('Time (years)','FontSize',16);
title('With Control','FontSize',16);
axis([0 inf 0 inf]);

set(gcf,'position',[442,634,865,333])
print(gcf,'ControlComparison.png','-dpng','-r300')
