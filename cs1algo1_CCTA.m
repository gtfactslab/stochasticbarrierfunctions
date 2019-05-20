% Example 1:    Verification Algorithm (i.e., Computes upper bound on
%               probability of failure.)
%
% Paper:    Verification  and  Control  for  Finite-Time  Safety  of Stochastic
%           Systems  via  Barrier  Functions (CCTA 2019)
%
% Author:   Cesar G. Santoyo
% Date:     May 17th, 2018
% 
% If you have any questions, e-mail the author at csantoyo@gatech.edu
%
% Description:  This script runs the verificcation algorithm for the 1-D
%               example of this paper. This paper requires the installation
%               of SOSTOOLS & SDPT3. This algorithm also has a dependency
%               with the file "example1algo.m". Generally speaking, this
%               algorithm performs a grid search over the variable alpha.

clc; close all;

alpha = [.1:.1:4];    % Range of alpha values for AB(x) =< -alpha*B(x) + Beta
sigma = [0:.1:2];     % Range of noise values

rho = 1;    % Value over which we want to stay below
T = 1;      % Time horizon (max time)
x0 = 0;     % Initial condition
deg = 16;   % Degree of Polynomials we are searching

% Pre-allocate optimal values
probvalue = ones(1,length(sigma));
gamprob = ones(1,length(sigma));
alpha_st = zeros(1,length(sigma));
beta_st = zeros(1,length(sigma));
gam_st = zeros(1,length(sigma));
u_st = zeros(1,length(sigma));

for ii = 1:length(sigma)        % Loop over noise level values
    for jj = 1:length(alpha)    % Loop over alpha values  
        
        [returnprob, retalpha, retbeta, Bx, retgam, prog] = cs1_verif(alpha(jj), sigma(ii), rho, T, x0, 0, deg);
        
        % Check if problem is feasible. If so, store values & compute
        % probablities of failure.
        if probvalue(ii) > returnprob && returnprob >= 0 && prog.solinfo.info.pinf == 0 && prog.solinfo.info.dinf == 0
           clc; % command window clear command used to prevent cluttering
           probvalue(ii) = returnprob
           alpha_st(ii) = retalpha
           beta_st(ii) = retbeta
           gam_st(ii) = retgam

           if beta_st(ii)/alpha_st(ii)<= 1
               gamprob(ii) = 1 - (1 - gam_st(ii))*exp(-beta_st(ii)*T)
           elseif beta_st(ii)/alpha_st(ii)>= 1
               gamprob(ii) = (gam_st(ii) + (beta_st(ii)/alpha_st(ii))*(exp(beta_st(ii)*T) - 1))/exp(beta_st(ii)*T)
           end
           Barriers(ii) = Bx;
           
        end
    end
end

%% Create Plots
figure;
plot(sigma,probvalue, 'LineWidth', 2, 'DisplayName', "$B(x) = x^{" + string(deg) + "}$")
hold on;
plot(sigma,gamprob, 'LineWidth',2, 'DisplayName','$\gamma$')
plot(sigma,count/draws, 'LineWidth', 2,'DisplayName',[string(draws) + " draw Monte Carlo"]); %% can comment out
grid on;
xlabel('g(x)','Interpreter','latex','Fontsize',14);
ylabel('Probabilities','Interpreter','latex','Fontsize',14);
title('Probabilities vs. g(x)','Interpreter','latex','Fontsize',16)
num = string(deg);
l = legend('show');
set(l, 'Interpreter','latex')
grid minor
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')
ylim([0 1])

figure;
plot(sigma, alpha_st, 'LineWidth', 2)
grid on;
xlabel('g(x)','Interpreter','latex','Fontsize',14);
ylabel('$\alpha^*$','Interpreter','latex','Fontsize',14)
title('$\alpha^*$ vs. g(x)','Interpreter','latex','Fontsize',16)
grid minor
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')

figure;
plot(sigma, beta_st, 'LineWidth', 2)
grid on;
xlabel('g(x)','Interpreter','latex','Fontsize',14);
ylabel('$\beta^*$','Interpreter','latex','Fontsize',14);
title('$\beta^*$ vs. g(x)','Interpreter','latex','Fontsize',16)
grid minor
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')

%% Plot level sets
% Description:  This section plots the level sets of the algorithm for
%               visual verfication. This portion of the code is optional
%               can be modified. This will help demonstrate the 
% x = [-2:.05:2];
% 
% Bvals = subs(Barriers);
% figure;
% plot(x, Bvals,'Linewidth',2)
% grid on
% grid minor
% xlabel('$x_1$','Interpreter', 'latex')
% ylabel('$B(x)$','Interpreter', 'latex')
% set(gcf,'color','w');
% set(gca,'TickLabelInterpreter','latex')
% 
% patch([-4 -1 -1 -4],[-3 -3 3 3], 'red', 'FaceAlpha', .2,'EdgeColor','black')
% patch(-1*[-4 -1 -1 -4],[-3 -3 3 3], 'red', 'FaceAlpha', .2,'EdgeColor','black')
% lim = 2.0;
% xlim([-lim lim])
% ylim([-lim lim])
% hold on;
% 
% % Plot initial set
% xc1 = 0;
% xc2 = 0;
% r = .2;
% hold on;
% th = 0:pi/50:2*pi;
% x1unit = r * cos(th) + xc1;
% x2unit = r * sin(th) + xc2;
% plot(x1unit,x2unit,'g--','LineWidth',1.5);
% plot(0,0,'r*','Linewidth',2)
% xlabel('$x_1$','Interpreter', 'latex')
% ylabel('$x_2$','Interpreter', 'latex')
% set(gcf,'color','w');
% set(gca,'TickLabelInterpreter','latex')
