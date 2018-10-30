clc; clear all; 

params.alpha = 216;%0.01;      
params.alpha0 = 0.216;%1;                     
params.n = 2;                        
params.beta = 0.3;

%% Initial state
tspan = [0, 1000]; %seconds
x0    = [0, 30, 35, 0, 0, 0];     %mRNA1,mRNA2,mRNA3, protein1, protein2, protein3

%% Specify reaction network
propensity_fcn = @propensities_2state;
% stoich_matrix = [ 1  0    %transcription
%                   0  1    %translation
%                  -1  0    %mRNA decay
%                   0 -1 ]; %protein decay
              
stoich_matrix = [ 1  0  0  0  0  0     %m1 transcription
                   0  1  0  0  0  0   %m2 transcription
                   0  0  1  0  0  0    %m3 transcription
                  -1  0  0  0  0  0    %m1 decay
                   0 -1  0  0  0  0      %m2 decay
                   0  0  -1 0  0  0    %m3 decay
                   0  0  0  1  0  0   %m1->p1 translation
                   0  0  0  0  1  0  %m2->p2 translation
                   0  0  0  0  0  1  %m3->p3 translation
                   0  0  0  -1 0  0      %p1 decay
                   0  0  0  0 -1  0      %p2 decay
                   0  0  0  0  0 -1];   %p3 decay



%% Initialize
MAX_OUTPUT_LENGTH = 1000000;
output_fcn = [];
%num_rxns = size(stoich_matrix, 1);
num_rxns = size(stoich_matrix, 1);
num_species = size(stoich_matrix, 2);
T = zeros(MAX_OUTPUT_LENGTH, 1);
X = zeros(MAX_OUTPUT_LENGTH, num_species);
T(1)     = tspan(1);
X(1,:)   = x0;
rxn_count = 1;

%% MAIN LOOP
while T(rxn_count) < tspan(2)        
    % Calculate reaction propensities
    a  = propensity_fcn(X(rxn_count,:), params);
    
    % Sample time-to-fire for each reaction channel
    % tau is the smallest time-to-fire and mu is the index of the channel
    r = rand(num_rxns,1);
    taus = -log(r)./a;
    [tau, mu] = min(taus);
    
%     tau = -log(rand(1))/a;
%     mus = randperm(num_rxns);
%     mu=  mus(1);

    
    % Update time and carry out reaction mu
    T(rxn_count+1)   = T(rxn_count)   + tau;
    X(rxn_count+1,:) = X(rxn_count,:) + stoich_matrix(mu,:);   
    rxn_count = rxn_count + 1; 
    

end  

% Return simulation time course
t = T(1:rxn_count);
x = X(1:rxn_count,:);
if t(end) > tspan(2)
    t(end) = tspan(2);
    x(end,:) = X(rxn_count-1,:);
end    


%% Plot time course
figure();
t_rescale= t*2.89;
X3proteins=[x(:,4).*40,x(:,5).*40,x(:,6).*40];
stairs(t_rescale,X3proteins); set(gca,'XLim',tspan);
xlabel('time (s)');
ylabel('molecules');

%%
function a = propensities_2state(x, params)
% Return reaction propensities given current state x
mRNA1    = x(1);
mRNA2    = x(2);
mRNA3    = x(3);
protein1 = x(4);
protein2 = x(5);
protein3 = x(6);

a = [params.alpha/(1+protein3^params.n)+params.alpha0;            %transcription
     params.alpha/(1+protein1^params.n)+params.alpha0;       %translation
     params.alpha/(1+protein2^params.n)+params.alpha0;       %mRNA decay
     mRNA1;
     mRNA2;
     mRNA3;
     params.beta*mRNA1;
     params.beta*mRNA2;
     params.beta*mRNA3;
     params.beta*protein1;
     params.beta*protein2;
     params.beta*protein3];  
end


