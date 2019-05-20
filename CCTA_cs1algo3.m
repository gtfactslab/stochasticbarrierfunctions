% Example 1:    Control Search Algorithm (i.e., Computes upper bound on
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
% Description:  This script performs a pseu-binary search for the the 1-D
%               controller of the form u(x) = -kx.

clc; clear all; 
degarr = [6 8 10 12 14 16];
% Add path to find other function files

for degcount = 1:length(degarr)
    alpha = [0:.1:5];    % Range of alpha values for AB(x) =< -alpha*B(x) + Beta
    gx = [1:.2:2];       % Range of noise values

    rho = 1;            % Value over which we want to stay below (i.e. safe set condition) 
    T = 1;              % Time horizon (max time)
    x0 = 0;             % Initial condition
    deg = degarr(degcount);            % Degree of Polynomials Barrier we are searching

    krange = [1:.1:50];                        % Range of input gains to search over
    offsetcount = 4;
    % Pre-allocate optimate values
    probvalue = ones(length(gx),floor(log(length(krange))) + offsetcount);             % Probability values for each noise lvl
    alpha_st = zeros(length(gx), floor(log(length(krange))) + offsetcount);             % Optimal alpha for each noise lvl
    beta_st = zeros(length(gx), floor(log(length(krange))) + offsetcount);              % Optimal beta for each noise lvl
    u_st = zeros(length(gx), floor(log(length(krange))) + offsetcount);                 % Optimal control input for each noise level
    % Initialize State variable x
    syms x

    prob_thresh = .30;                          % Desired probability of failure for noise level


    k_st = -999*ones(length(gx), floor(log(length(krange))) + offsetcount);

    for ii = 1:length(gx)   % Loop through noise levels
        % Initialize control variable & auxiliary variables
        infcount = 1;
        lo = 1;                                     % Lower bound of binary search  (index value)                       
        up = length(krange);                        % Upper bound of binary search  (index value)

        % loop while the probability value difference is not within some range
        while ~(probvalue(ii) - prob_thresh > -1e-2 && probvalue(ii) - prob_thresh < 0)
            % If statement to catch infinte loops. Iteration count limited to
            % log(n) where n = elements in gain grid range
            if infcount > floor(log(length(krange))) + offsetcount
                if m > 0
                    k_st(ii, end) = -krange(m);
                else
                    k_st(ii, end) = 0;
                end

                break;  % End while loop because iterating for too long
            end

    %         probvalue(ii) = 1;          % Re-assign the probability to 1
            m = floor((lo + up)/2);     % Assign mid way value

            % Assign control input based on control gain ranges
            if m > 0
                u = -krange(m)*x;
            else
                u = 0;
            end

            % Grid search through potential alpha candidates
            for jj = 1:length(alpha)
                [returnprob, retalpha, retbeta, Bx] = cs1_verif(alpha(jj), gx(ii), rho, T, x0, u, deg)
                if probvalue(ii, infcount) > returnprob
                   clc
                   Barrfunc(ii, infcount) = Bx;
                   probvalue(ii, infcount) = returnprob
                   alpha_st(ii, infcount) = retalpha
                   beta_st(ii, infcount) = retbeta

                   if m > 0
                        k_st(ii, infcount) = -krange(m)
                   else
                        k_st(ii, infcount) = 0;
                   end
                end
            end

            % Check conditions for binary search to assign new lower & upper
            % bounds
            if probvalue(ii, infcount) > prob_thresh          
                lo = m + 1;
            elseif probvalue(ii, infcount) < prob_thresh     
                up = m - 1;
            end

            % Check to see if while-loop is about to exit. If so, save values
            if (probvalue(ii, infcount) - prob_thresh > -1e-2 && probvalue(ii, infcount) - prob_thresh < 0) 
                if m > 0
                    k_st(ii, infcount) = -krange(m);
                else
                    k_st(ii, infcount) = 0;
                end

                break;
            end

            % If prob value is very close to zero, break while loop
            if probvalue(ii, infcount) < .01 && abs(k_st(ii, infcount))< .05
                break;
            end

            infcount = infcount + 1;
        end
    end


    % Plot figures of probability bound
    probplot = zeros(1,length(gx));
    alphaplot = zeros(1,length(gx));
    betaplot = zeros(1,length(gx));
    kplot = zeros(1,length(gx));
    for ii = 1:length(gx)
        for jj = 1:size(probvalue,2)
            if probvalue(ii, jj) < prob_thresh && probvalue(ii, jj) > probplot(ii)
                probplot(ii) = probvalue(ii, jj);
                alphaplot(ii) = alpha_st(ii,jj);
                betaplot(ii) = beta_st(ii,jj);
                kplot(ii) = k_st(ii,jj);
            end
        end
    end
    save([string(datetime('now','format','DDMMMyy_HHmmss'))+'_controlsearch_deg' + string(degarr(degcount))]);
end
%% Plot Figures
figure;
plot(gx,probplot, 'LineWidth', 2)
grid on;
xlabel('g(x)','Interpreter','latex','Fontsize',14);
ylabel('Probability','Interpreter','latex','Fontsize',14);
title('Probabilities vs. g(x)','Interpreter','latex','Fontsize',16)
legend(['B(x) = x^{' + string(deg) + '}'])
ylim([0 1]);
grid minor
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')

figure;
plot(gx, alphaplot, 'LineWidth', 2)
grid on;
xlabel('$\sigma (x)$','Interpreter','latex','Fontsize',14);
ylabel('$\alpha^*$','Interpreter','latex','Fontsize',14)
title('$\alpha^*$ vs. $\sigma (x)$','Interpreter','latex','Fontsize',16)
grid minor
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')

figure;
plot(gx, betaplot, 'LineWidth', 2)
grid on;
xlabel('$\sigma (x)$','Interpreter','latex','Fontsize',14);
ylabel('$\beta^*$','Interpreter','latex','Fontsize',14);
title('$\beta^*$ vs. $\sigma (x)$','Interpreter','latex','Fontsize',16)
grid minor
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')

figure;
plot(gx, -1*kplot, 'LineWidth', 2)
grid on;
xlabel('$\sigma$','Interpreter','latex','Fontsize',14);
ylabel('$k^*$','Interpreter','latex','Fontsize',14);
% title('Control Gain vs. $\sigma (x)$','Interpreter','latex','Fontsize',16)
grid minor
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')
hold on;
% figure;
% dt = delaunayTriangulation(gx', -kplot') ;
% tri = dt.ConnectivityList ;
% xi = dt.Points(:,1) ; 
% yi = dt.Points(:,2) ; 
% F = scatteredInterpolant(gx', -kplot', probplot');
% zi = F(xi,yi) ;
% trisurf(tri,xi,yi,zi) 
% view(2)
% shading interp
% colormap winter
% h = colorbar;
% ylabel(h, 'Probability','Interpreter','latex','Fontsize',14)
% xlabel('$g(x)$','Interpreter','latex','Fontsize',14);
% ylabel('$k^*$','Interpreter','latex','Fontsize',14);
% zlabel('$Probability$','Interpreter','latex','Fontsize',14);

