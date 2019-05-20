% Example 1:    Monte Carlo Simulation of 1-D Dynamics using SDE & 
%               
% Paper:    Verification  and  Control  for  Finite-Time  Safety  of Stochastic
%           Systems  via  Barrier  Functions (CCTA 2019)
%
% Author:   Cesar G. Santoyo
% Date:     May 17th, 2018
% 
% If you have any questions, e-mail the author at csantoyo@gatech.edu
%
% Description:  This script runs & plots the Monte Carlo simulation. Two
%               techniques are used to plot the monte carlo dynamics. There
%               is a discrete time formulation of the standard Weiner
%               process and the continuous time version which uses the SDE
%               toolbox from MATLAB. 

%% DT Modeling of Continuous Time Dynamics
clear all; clc;
sigma = [0:.1:2];   % Noise Level
mean = 0;           % Mean of our normal distribution
stddev = 1;         % Standard deviation of our normal distribution
draws = 5000;       % Number of Monte Carlo draws (design parameter)
dt = .01;           % Time Step

%Pre-allocate memory
noise = zeros(length([0:dt:1]),draws,length(sigma));
count = zeros(1,length(sigma));

for k = 1:length(sigma)
    for  i = 1:draws
        for j = 1:length([0:dt:1])
            if j == 1
                x(j) = 0;   % Initial Condition
            else
                noise(j,i,k) = randn(1);    % Draw Noise
                x(j) = (1 - dt).*x(j - 1) + sqrt(dt)*sigma(k)*noise(j,i,k); % DT Modeling of Standard Weiner Process
            end
        end
        if ~isempty(x(x>=1))
            count(k) = count(k) + 1;
        end
    end
end

figure;
plot(sigma,count/draws, 'LineWidth', 2,'DisplayName',[string(draws) + " draw Monte Carlo"]);
grid on;
grid minor
xlabel('$\sigma$', 'Interpreter', 'latex');
ylabel('Probability', 'Interpreter', 'latex');
ylim([0 1]);
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex','FontSize', 14)
l = legend('show');
set(l,'Interpreter','latex');

%% Monte Carlo Simulation using SDE Toolbox 

% NOTE: This Monte Carlo requires access to MATLAB's SDE Toolbox

% dXt=F(t,Xt)dt+G(t,Xt)dWt

%F is an NVARS-by-1 vector-valued drift-rate function
%G is an NVARS-by-NBROWNS matrix-valued diffusion-rate function.
clc; clear all;

sigma = [0:.1:2];
dt = .01;
nPeriods = 100;
draws = 5000;
quantity_true = zeros(draws,length(sigma));
count = zeros(1, length(sigma));
for sigma_idx = 1:length(sigma)
    for draws_idx = 1:draws
        driftrate = drift(0, -1);
        diffrate = diffusion(0,sigma(sigma_idx));

        SDE = sde(driftrate, diffrate,'StartState',0);

        [Paths,Times,Z] = simulate(SDE, nPeriods, 'DeltaTime' , dt);
        if numel(Paths(Paths>=1)) > 0
            count(sigma_idx) = count(sigma_idx) + 1;
        end
    end
end

prob_count = count/draws;
figure;
plot(sigma,prob_count,'LineWidth',2,'DisplayName', 'SDE Sim (CT)');
grid on;
xlabel('g(x)');
ylabel('Probability');
ylim([0 1]);
