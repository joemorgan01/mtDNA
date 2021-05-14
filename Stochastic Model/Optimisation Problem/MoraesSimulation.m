function [SSE2]=MoraesSimulation(Parameters)
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

Error = NaN(4,1);
Error(1) = MoraesRedLine(3,2) - normalisedcopyno(1);
Error(2) = MoraesRedLine(4,2) - normalisedcopyno(7*24);
Error(3) = MoraesRedLine(5,2) - normalisedcopyno(15*24);
Error(4) = MoraesRedLine(6,2) - normalisedcopyno(30*24);

SqError=Error.^2;
SSE2=sum(SqError,'all');
end

