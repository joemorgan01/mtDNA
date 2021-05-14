function [SSE2]=MoraesFinal(Parameters)
%% Define units of time and set simulation end time
year=365*24*3600;
day=24*3600;
rxnT=30*day;

%% Set number of Simulations to perform
RunNo=1000;

%% Sample Time 
sampletime=linspace(0,rxnT,24*365*rxnT/year);

%% Initialise Array Sizes
Yplot200=NaN(2, length(sampletime));
tplot=NaN;
Results=NaN(length(sampletime),2 , RunNo);
TResults=NaN(1,length(sampletime), RunNo);
COXDefs=NaN(length(sampletime),1,RunNo);
Inlet=NaN(1,2,RunNo);
Rates=NaN(1,5,RunNo);
Mol=NaN(1,1,RunNo);

%% Determine simulation condition by sampling from probability distribution
% funcitons
for i=1:RunNo
    MoleculeNo=MoraesCopyNo(Parameters);
    Y=InitialConditionsv2(MoleculeNo);
    C0=Rates4Optimisation(Parameters);
    Rates(:,:,i)=C0;
    Inlet(:,:,i)=Y;
    Mol(:,:,i)=MoleculeNo;
end

%% Run Simulation
for i=1:RunNo
    InitialCondition=Inlet(:,:,i);
    Rate=Rates(:,:,i);
    [tplot,Yplot200,COXDef] = mtDNAModelMoraes([InitialCondition(1) InitialCondition(2)], sampletime, Rate, Parameters);
    Results(:,:,i) = Yplot200;
    TResults(:,:,i)=tplot;
    COXDefs(:,:,i)=COXDef;
    i=i;
end


load('MoraesBlueLine');
load('MoraesRedLine');

copynumber=(Results(:,2,:)+Results(:,1,:));
avcopyno=mean(copynumber,3);
normalisedcopyno=avcopyno./Parameters(2);
%figure(1)
 subplot(1,2,2)
avcopyno=mean(copynumber,3);
copyno95=prctile(copynumber,95,3);
copyno5=prctile(copynumber,5,3);
stairs((sampletime./day)+15,100*avcopyno./Parameters(2),'k-','LineWidth',3);
hold on
stairs((sampletime./day)+15,100*copyno95./Parameters(2),'k--','LineWidth',0.5);
hold on
stairs((sampletime./day)+15,100*copyno5./Parameters(2),'k--','LineWidth',0.5);
hold on
errorbar(MoraesRedLine(:,1),100*MoraesRedLine(:,2),100.*[0 0.0312 0.0312 0.0424 0.0981 0.1],'ko','LineWidth',1,'MarkerFaceColor','k')


grid on
ylabel('Copy Number as % of Setpoint','FontSize',16);
xlabel('Time (days)','FontSize',16);
axis([0 inf 0 inf]);

Error = NaN(4,1);
Error(1) = MoraesRedLine(3,2) - normalisedcopyno(1);
Error(2) = MoraesRedLine(4,2) - normalisedcopyno(7*24);
Error(3) = MoraesRedLine(5,2) - normalisedcopyno(15*24);
Error(4) = MoraesRedLine(6,2) - normalisedcopyno(30*24);

SqError=Error.^2;
SSE2=sum(SqError,'all');
end

