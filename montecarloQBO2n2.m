% TEAM ASSIGNMENT 1: MONTECARLO METHOD
function [end_t, n, G]=montecarlo(T,P,B,maxt)

%% Properties
species ={'H2','O2','H2O','H','O','OH'};

n1 = [2;1];
n2 = [0;0;0;0]; % Initial number of moles

A1 = [2 0 ; ...
      0 2 ]; % ≠ types and number of atoms in each specie 
A2 = [2 1 0 1; ...
      1 0 1 1]; % ≠ types and number of atoms in each specie 

Gold = 0;% Gibbs energy. Needs to be <= 0.
tic % Starting time for calculation 

% for i = 1:1e5
current_t = 0;
iter = 0;

while current_t < maxt-0.001
    
    n2old = n2; % save previous values 
    n1old = n1; %  "      "       "

    n2 = max(B)*rand(length(n2),1); % rand is always 0 <=x <= 1. So n2 + rand > 0
    n1 = A1 \ (B - A2*n2); % n1 will be <0 if we consume too much mol
    n = [n1; n2]; % mol by species to satisfy stochiometric coeff

    current_t=toc;
    iter=iter+1;
    
    if min(n) >= 0 % the number of atoms cannot be negative
        [~,~,~,~,~,~,~,G,~]=hgsprop(species,n,T,P);  
            if  Gold>G       
                Gold=G; % We keep the lowest calculated value of G and we search for an even lower value         
            end
    else
        n = [n1old;n2old];
        continue
    end
    
end

iter
end_t=toc; % Ending time of calculation