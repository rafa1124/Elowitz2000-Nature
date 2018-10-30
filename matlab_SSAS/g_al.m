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
                   0  0  0  1  0  0  %m1->p1 translation
                   0  0  0  0  1  0%m2->p2 translation
                   0  0  0  0  0  1 %m3->p3 translation
                   0  0  0  -1 0  0      %p1 decay
                   0  0  0  0 -1  0      %p2 decay
                   0  0  0  0  0 -1];   %p3 decay



%% Initialize
MAX_OUTPUT_LENGTH = 1000000;
output_fcn = [];
%num_rxns = size(stoich_matrix, 1);
num_species = size(stoich_matrix, 2);
T = zeros(MAX_OUTPUT_LENGTH, 1);
X = zeros(MAX_OUTPUT_LENGTH, num_species);
T(1)     = tspan(1);
X(1,:)   = x0;
rxn_count = 1;

%% MAIN LOOP
while T(rxn_count) < tspan(2)        
    % Calculate reaction propensities
    a = propensity_fcn(X(rxn_count,:), params);
    
    % Sample earliest time-to-fire (tau)
    a0 = sum(a);
    r = rand(1,2);
    tau = -log(r(1))/a0; %(1/a0)*log(1/r(1));
    
    % Sample identity of earliest reaction channel to fire (mu)
    [~, mu] = histc(r(2)*a0, [0;cumsum(a(:))]); 
    
    % ...alternatively...
    %mu = find((cumsum(a) >= r(2)*a0), 1,'first');
    
    % ...or...
    %mu=1; s=a(1); r0=r(2)*a0;
    %while s < r0
    %   mu = mu + 1;
    %   s = s + a(mu);
    %end

    if rxn_count + 1 > MAX_OUTPUT_LENGTH
        t = T(1:rxn_count);
        x = X(1:rxn_count,:);
        warning('SSA:ExceededCapacity',...
                'Number of reaction events exceeded the number pre-allocated. Simulation terminated prematurely.');
        return;
    end
    
    % Update time and carry out reaction mu
    T(rxn_count+1)   = T(rxn_count)   + tau;
    X(rxn_count+1,:) = X(rxn_count,:) + stoich_matrix(mu,:);    
    rxn_count = rxn_count + 1;
    
    if ~isempty(output_fcn)
        stop_signal = feval(output_fcn, T(rxn_count), X(rxn_count,:)');
        if stop_signal
            t = T(1:rxn_count);
            x = X(1:rxn_count,:);
            warning('SSA:TerminalEvent',...
                    'Simulation was terminated by OutputFcn.');
            return;
        end 
    end
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

