% Example 2:    Monte Carlo Simulation of 2-D Dynamics using SDE & 
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
%               process for Case Study 2 in this paper. 

%% Euler Method Monte Carlo Stochastic van der Pol 
clear all; close all;
T = 2.0;                       % Time Horizon
dt = .01;                      % Time step
iterations = floor(T/dt);      % Number of iterations needed 
mcidx = 5000;                  % Max number of MonteCarlo draws

muval = 1;              % van der Pol damping value
meanval = 0;            % Noise distribution mean
stddev = 1;             % Standard deviation of noise distribution
gx = [.5:.1:1.5];          % Array of Noise levels

% Pre-allocate array to count number of iterations
unsafecount = zeros(1, length(gx));
ux = [];
for sigidx = 1:length(gx)       % Loop through noise levels
    % Pre-allocate state and state time derivative
    x = zeros(2,iterations);
    dx = zeros(2,iterations);
    
    for ii = 1:mcidx            % Loop through Monte Carlo iterations
        
        x(:,1) = [-2;0];        % Initial Condition

            % Stable at the origin dynamics
        dx(1,1) = x(2,1)*dt;            
        x1 = x(1,1);
        x2 = x(2,1);
        ux(1) = (87058438113*x1)/137438953472 - (175733967931*x2)/274877906944 + x1*((275195185673*x2)/549755813888 - (351576535553*x1)/549755813888 + 174117523085/274877906944) + x2*((17199844495*x1)/34359738368 + (4590613*x2)/549755813888 - 351465519453/549755813888) - 351271316573/549755813888;
        dx(2,1) = (- x(1,1) - x(2,1) - .5*x(1,1)^3 + ux(1))*dt;


        for jj = 2:iterations   % Loop through time steps for complete time horizon
            x(1,jj) = x(1,jj - 1) + dx(1,jj - 1);
            x(2,jj) = x(2,jj - 1) + dx(2,jj - 1);
            
            x1 = x(1,jj);
            x2 = x(2,jj);
            
            % One of the controllers synthesized using SOS
%             ux(jj) = (87058438113*x1)/137438953472 - (175733967931*x2)/274877906944 + x1*((275195185673*x2)/549755813888 - (351576535553*x1)/549755813888 + 174117523085/274877906944) + x2*((17199844495*x1)/34359738368 + (4590613*x2)/549755813888 - 351465519453/549755813888) - 351271316573/549755813888;
            ux(jj) = 0;
            
            % Stable at the origin dynamics
              dx(1,jj) = x(2,jj)*dt;
              dx(2,jj) = (- x(1,jj) - x(2,jj) -.5*x(1,jj)^3 + ux(jj))*dt + sqrt(dt)*gx(sigidx)*normrnd(meanval,stddev);

        end

        
        if sum(x(2,:) >= 2.25) > 0
            unsafecount(sigidx) = unsafecount(sigidx) + 1;
        end
    end
end
% Plot solution
figure;
plot(gx, unsafecount/mcidx,'k','LineWidth', 2)
xlabel('$\sigma$','Interpreter', 'latex')
ylabel('Probability','Interpreter', 'latex')
title(['van der Pol Oscillator Probabilty vs. Noise Level (T = '+string(T)+ ' sec)'],'Interpreter', 'latex')
grid on;
grid minor
ylim([0 1]);
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')
figtitle = strcat('montecarlo',string(mcidx),'draws');
savefig(char(figtitle));
print('montecarlo','-dpng','-r300')
% figure;
% plot(x(1,1:jj-1), x(2,1:jj-1),'k','LineWidth', 2)
% xlabel('$x_1$','Interpreter', 'latex')
% ylabel('$x_2$','Interpreter', 'latex')
% title(['van der Pol Oscillator, $\mu = ' + string(muval)+'$'+' (T = ' + string(T) + ' sec)'],'Interpreter', 'latex')
% grid on;
% grid minor
% set(gcf,'color','w');
% set(gca,'TickLabelInterpreter','latex')