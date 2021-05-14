function [SSE1,ks2stat]=RochaFinal(Parameters)
%% Define units of time and set simulation end time
year=365*24*3600;
day=24*3600;
rxnT=80*year;


%% Set number of Simulations to perform
RunNo=100;

%% Sample Time 
sampletime=linspace(0,rxnT,12*rxnT/year);

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
    MoleculeNo=RochaCopyNo(Parameters);
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
    [tplot,Yplot200,COXDef] = mtDNAModelRocha([InitialCondition(1) InitialCondition(2)], sampletime, Rate, Parameters);
    Results(:,:,i) = Yplot200;
    TResults(:,:,i)=tplot;
    COXDefs(:,:,i)=COXDef;
    i=i;
end


load('UnscaledRochasCN');
% Need to scale data here to have median value of Parameters(1)
ScaledMedian=Parameters(1);
UnscaledMedian=median(UnscaledRochasCN,'all');
ScaleFactor=ScaledMedian/UnscaledMedian;
ScaledRochas=UnscaledRochasCN.*ScaleFactor;
Med=median(ScaledRochas,'all')
copynumber=(Results(:,2,:)+Results(:,1,:));
CopyNumber=NaN(length(sampletime)*RunNo,1);
CopyNumber(1:length(sampletime))=copynumber(:,:,1);

for i=1:(RunNo-1)
    CopyNumber((length(sampletime)+1)*i:(((length(sampletime)+1)*i)+(length(sampletime)-1)))=copynumber(:,:,i+1);
end

[SimDensity, x1]=ksdensity(CopyNumber);
[RochaDensity,x2]=ksdensity(ScaledRochas,'BandWidth',25);
Error=SimDensity-RochaDensity;
SqError=Error.^2;
SSE1=sum(SqError,'all');

[~,~,ks2stat]=kstest2(SimDensity,RochaDensity);

histogram(CopyNumber,20,'Normalization','probability','EdgeColor','none','FaceColor',[0.5 0.5 0.5]);
hold on
histogram(ScaledRochas,100,'Normalization','probability','EdgeColor','none','FaceColor',[0 0 0]);
xlabel('Copy Number','FontSize',16);
ylabel('Probability','FontSize',16);
axis([0 2000 0 inf]);
legend('Model','Rocha et Al','FontSize',14);
grid on
end

