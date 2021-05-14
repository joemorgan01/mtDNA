clear all;
close all;
% Define Time and Reaction Duration
year=365*24*3600;
day=24*3600;
rxnT=80*year;

% Specify Number of Simulations to Perform
RunNo=1000;

% Set sample time
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
    [tplot,Yplot200,COXDef] = mtDNAModelv2wControl([InitialCondition(1) InitialCondition(2)], sampletime, Rate);
    Results(:,:,i) = Yplot200;
    TResults(:,:,i)=tplot;
    COXDefs(:,:,i)=COXDef;
    i=i
end
disp('Done');
