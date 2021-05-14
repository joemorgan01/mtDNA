function [dCdt] = InheritedODEModel (~,~)
%----------------------------
% Equation          | Rate
%-------------------|--------
% A --> 2A          |   C1
% B --> 2B          |   C2
% A --> 0           |   C3
% B --> 0           |   C4
% A --> A + B       |   C5


%       A  B
Stoich=[1  0    % R1
        0  1    % R2
        -1 0    % R3
        0 -1    % R4
        0  1];  % R5
    
ReactionRate = [8.0225e-7,... % C1
                8.0225e-7,... % C2 
                8.0225e-7,... % C3
                8.0225e-7,... % C4 
                1.157e-10];   % C5
            
dCdt=Stoich'*ReactionRate';
end