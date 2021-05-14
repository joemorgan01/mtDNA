function [tplot, Yplot, COXDef] = mtDNAModelMoraes ( Y0, timeSamples, C0, Parameters)
% function to sumulate the network of reactions which make up the
% synthesis, degredation, and mutation of mtDNA
% The inputs of this function are the initial number of molecules and the desired 
% sampling times.
%% Key
% Wild Type = A
% Mutant    = B
%--------------------------------------------------------------------------
%Y0           | number of molecules of each species at time t=0
%             |
%timeSamples  | vector of times at which we would like to know the number
%             | of each species
%--------------------------------------------------------------------------
%M = 5         Number of reaction pathways
%N = 2         Number of molecular species considered

%% Step 0
% Initialise number of modelcules of species A and B. In general, a vector
% of length N of initial population numbers
%----------------------------
% Equation          | Rate
%-------------------|--------
% A --> 2A          |   c1
% B --> 2B          |   c2
% A --> 0           |   c3
% B --> 0           |   c4
% A --> A + B       |   c5

%% Time Sampling Code
day=24*3600;
year=365*day;
obs_index = 1;
n_times_to_observe = length ( timeSamples ); % number of sample times

Y=Y0;

Yplot = NaN( n_times_to_observe, length (Y));
COXDef= NaN(n_times_to_observe,1);
error=NaN(n_times_to_observe,1);
t=0;
e=0;
tplot = NaN;

%% Step 1
while obs_index<= n_times_to_observe
    
    if (t>=timeSamples(obs_index))
        Yplot(obs_index,:) = Y;
        tplot(obs_index) = t;
        error(obs_index,:)= e;
        
        % If at a sample time the mutation % is above the critical
        % threshold (80%), then record the COX deficiency as '1'
        % Otherwise, record COX deficiency as '0'.
        if (Y(2)/sum(Y)) > 0.8
            COXDef(obs_index,:)=1;
        else
            COXDef(obs_index,:)=0;
        end
        
        obs_index = obs_index + 1;
        continue;
    end
    
% Controller
    Setpoint=Parameters(2);  
    e=Setpoint-sum(Y);
    if e<0 %|| abs(e)<(0.1*Setpoint)
        MV=0;
    else
        Kc=Parameters(4);
        MV=e*Kc;
    end

%     MV=0;
    
 % stochastic reaciton constants
 % c1-4 = ln(2)/half life (10 days)
 
c=[ (C0(1)+MV),... 
    (C0(2)+MV),...
    C0(3),...
    C0(4),...
    C0(5)]; % seconds
    % define the hazard of each event happening, since all of these
    % reactions are first-order, this is the same as the rate equation
    h=[Y(1),...
       Y(2),...
       Y(1),...
       Y(2),...
       Y(1)];
   a=h.*c;
   a0= sum(a);
   
   % Probability of each event happening
   
   p=a./a0;
   p=cumsum(p);
   
   %% Step 2
% Determine next sample time using exponential distrubution function with an average
% value of 1/ (total hazard)
   tau=exprnd(1/a0);
   t=t+tau;
% Which reaction occurs is determined by a uniformly distributed random
% number within the range 0-1
   mu = rand;
   
   %% Step 3

% Adjust population of species depending on which reaction is selected by
% mu
   if (0<=mu)&&(mu<=p(1))
       Y(1) = Y(1) + 1;
   end
   if (p(1)<mu)&&(mu<=p(2))
      Y(2) = Y(2) + 1;
   end
   if (p(2)<mu)&&(mu<=p(3))
           Y(1) = Y(1) - 1;
   end
   if (p(3)<mu)&&(mu<p(4))
           Y(2) = Y(2) - 1;
   end
   if (p(4)<mu)&&(mu<=p(5))
           Y(2) = Y(2) + 1;
   end

end
      
end
        


