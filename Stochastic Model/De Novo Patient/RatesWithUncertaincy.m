function c = RatesWithUncertaincy()
day=24*3600;
% Determine the half life of the mitochondrial DNA by sampling from a
% gaussian distribution with a mean value of 260 days.

r=normrnd(260,1);


% convert half life in days to reaciton rate in seconds

rate=log(2)/r; % days
rate=rate./(day); %seconds

% Array of reaction rates to be used in the main mtDNA model

c=[rate,... % Synthesis of A rate
    rate,...% Synthesis of B rate
    rate,...% Degredation of A rate
    rate,...% Degredation of B rate
    1.157e-12]; % Rate of mutation (uncertanty not considered);
end