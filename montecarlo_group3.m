% TEAM ASSIGNMENT 1: MONTECARLO METHOD group 3
%Pol Gil, Núria Tirado, Tomás Deus and Quentin Boulogne
function [n]=montecarlofinalV(T,P,B,maxt)
%addpath('...\..\HGS_2020\')% Changing the path to your path to HGS_2020
tic
% Properties
species ={'H2','O2','H2O','H','O','OH'};
A1 = [2 0 ; ...
      0 2 ]; % â‰  types and number of atoms in each specie 
A2 = [2 1 0 1; ...
      1 0 1 1]; % â‰  types and number of atoms in each specie 
Gold = 0;% Gibbs energy. Needs to be <= 0.
current_t=0;
iter=0;
while current_t<maxt-0.001 % to do calculation for 10 sec (-0.001 to be below the 10 sec at the end
        n2 = max(B)*rand(4,1); % rand is always 0 <=x <= 1. So n2 => 0
        n1 = A1 \ (B - A2*n2); % n1 may be <0
        nt = [n1; n2]; % mol by species to satisfy stochiometric coeff
        if min(n1) >= 0 % the number of mol cannot be negative (n1 is the only values that could be <0)
            [~,~,~,~,~,~,~,G,~]=hgsprop(species,nt,T,P);% Computing G with hgprop
            if  Gold>G  % compute if G is lower than the lowest G computed on previous iterations
                n=nt;%return n for the lowest G value 
                Gold=G;% We keep the lowest calculated value of G and we search for an even lower value         
            end
        end
        current_t=toc;  % Current time
        iter=iter+1;    % Number of iteration 
end 
end_t=toc; % Ending time of calculation will be almost 10 s