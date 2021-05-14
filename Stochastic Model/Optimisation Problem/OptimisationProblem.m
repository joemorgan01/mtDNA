clear all
close all

% %Initial values for optimiser
% % Setpoint for Rocha Data
% Setpoint1=200;
% % Setpoint for Moraes Data
% Setpoint2=200;
% ReplicationRate=3.056e-8;
% ControllerGain=9e-9;
% 
% %% Initial Guess
% Parameters=[Setpoint1, Setpoint2, ReplicationRate, ControllerGain];
% 
% %% Define equality constraints
% Aeq =[];
% Beq=[];
% 
% %% Set upper and lower bound 
% LB =[199 0 2e-8 0];
% UB =[1e5 500 8e-6 8e-6];
% 
% %% Set inequality constraints 
% A =[0 0 -1 1];
% B = 0;
% 
% %% Call the optimiser
% [Parameters, TotalError] = fmincon(@(Parameters)OptimisationProblems(Parameters),...
%                                Parameters,A,B,Aeq,Beq,LB,UB);

load('OptimisedParameters');
%% Plot Results using the optimised parameters                           
[~]=Results(Parameters);
set(gcf,'position',[442,634,865,333])
print(gcf,'OptimisationFig.png','-dpng','-r300')
%% Functions
function [TotalError] = OptimisationProblems(Parameters)
%% Rocha et al
[SSE1,ks2stat]=RochaSimulation(Parameters);
%% Moraes
[SSE2]=MoraesSimulation(Parameters);
%% Error to minimise
TotalError=(SSE2+ks2stat)
end

function [TotalError] = Results(Parameters)
%% Rocha et al
[SSE1,ks2stat]=RochaFinal(Parameters);
%% Moraes
[SSE2]=MoraesFinal(Parameters);
%% Error to minimise
TotalError=SSE2+ks2stat;
end