function c = Rates4Optimisation(Parameters)

rate=normrnd(Parameters(3),0.15*Parameters(3));

% Array of reaction rates to be used in the main mtDNA model

c=[rate,... % Synthesis of A rate
    rate,...% Synthesis of B rate
    rate,...% Degredation of A rate
    rate,...% Degredation of B rate
    1.157e-10./50]; % Rate of mutation (uncertanty not considered);
end