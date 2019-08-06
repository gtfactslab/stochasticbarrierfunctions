% Main Files to run examples
% This file computes the probability bound based on this new form of
% infinitesimal generator
clc; close all;
syms x1 x2
alpha = [.1:.1:2];  % Range of alpha values for AB(x) =< -alpha*B(x) + Beta
gx = [.5:.1:1.5];     % Range of noise values
barrdeg =[14];
rho = 1;        % Value over which we want to stay below
T = 2.0;        % Time horizon (max time)
x0 = [-2 0];     % Initial condition
u = 0;        % Control input


for kk = 1:length(barrdeg)
    deg = barrdeg(kk);       % Degree of Polynomials we are searching
    
    % Pre-allocate optimate values
    probvalue = ones(1, length(gx));
    alpha_st = zeros(1, length(gx));
    beta_st = zeros(1, length(gx));
    gam_st = zeros(1,length(gx));
    gamprob = zeros(1,length(gx));

    for ii = 1:length(gx)
        for jj = 1:length(alpha)
                [returnprob, retalpha, retbeta, retgam, Bx, prog] = cs2_verif(alpha(jj), gx(ii), rho, T, x0, u, deg)
                if probvalue(ii) > returnprob && returnprob >= 0 && prog.solinfo.info.pinf == 0 && prog.solinfo.info.dinf == 0 
                   clc;
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
                elseif prog.solinfo.info.pinf == 1 || prog.solinfo.info.dinf == 1 
                    disp('Primal or dual infeasible problem.');
                end
        end
    end
    
save(['kushnernoctrl_polydeg' + string(deg) + '_' +  string(datetime('now','format','DDMMMyy_HHmmss'))])
end
%%
figure;
plot(gx,probvalue, 'LineWidth', 2)
hold on;
plot(gx, unsafecount/mcidx,'k','LineWidth', 2)
grid on;
xlabel('$\sigma$','Interpreter','latex','Fontsize',14);
ylabel('Probabilities','Interpreter','latex','Fontsize',14);
title(['Probabilities vs. $\sigma$ (T = '+ string(T) + ' sec)'],'Interpreter','latex','Fontsize',16)
legend(['B(x) = x^{' + string(deg) + '}'])
ylim([0 1]);
grid minor
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')

figure;
plot(gx, alpha_st, 'LineWidth', 2)
grid on;
xlabel('$\sigma$','Interpreter','latex','Fontsize',14);
ylabel('$\alpha^*$','Interpreter','latex','Fontsize',14)
title('$\alpha^*$ vs. $\sigma$','Interpreter','latex','Fontsize',16)
grid minor
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')

figure;
plot(gx, beta_st, 'LineWidth', 2)
grid on;
xlabel('$\sigma$','Interpreter','latex','Fontsize',14);
ylabel('$\beta^*$','Interpreter','latex','Fontsize',14);
title('$\beta^*$ vs. $\sigma$','Interpreter','latex','Fontsize',16)
grid minor
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')

%% Plot level sets
x1val = [-4:.05:4];
x2val = [-4:.05:4];

[X1 X2] = meshgrid(x1val, x2val);
x1 = X1;
x2 = X2;
Bvals = subs(Barriers(1));
% figure;
[cm,c]=contour(X1, X2, Bvals, [0 1 3 20],'b-.','LineWidth',1.2);
clabel(cm, c,'Color','k','FontSize', 14, 'Interpreter','latex','LabelSpacing', 210);
grid on
grid minor
xlabel('$x_1$','Interpreter', 'latex')
ylabel('$x_2$','Interpreter', 'latex')
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')

patch([-4 4 4 -4],[2.25 2.25 4 4], 'red', 'FaceAlpha', .2,'EdgeColor','black')
lim = 4.0;
xlim([-lim lim])
ylim([-lim lim])
hold on;

% Plot initial set
xc1 = -2;
xc2 = 0;
r = .1;
hold on;
th = 0:pi/50:2*pi;
x1unit = r * cos(th) + xc1;
x2unit = r * sin(th) + xc2;
plot(x1unit, x2unit,'g--', 'LineWidth', 2);
% plot(0,0,'r*','Linewidth',2)
xlabel('$x_1$','Interpreter', 'latex')
ylabel('$x_2$','Interpreter', 'latex')
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')


%% Plot Bound dependent on gamma
[value, idx] = min(probvalue);
for ii = 1:size(gam_st,2)
    gamvalbound(ii) = gam_st(1,ii);
    if 1 >=  beta_st(1, ii)/alpha_st(1, ii)
        probgam(ii) = 1 - ( 1 - gamvalbound(ii))*exp(-beta_st(1, ii)*T) ;
    elseif 1 <=  beta_st(1, ii)/alpha_st(1, ii)
        probgam(ii) =  (gamvalbound(ii) + (exp(beta_st(1, ii)*T) - 1)*(beta_st(1, ii)/alpha_st(1, ii)) )/(rho*exp(beta_st(1, ii)*T));
    else 
        probgam(ii) = (gamvalbound(ii + beta_st(1, ii)*T)/1);
    end
end

figure;
plot(gx, probgam);
grid on;
grid minor;
xlabel('Probability','Interpreter','latex');
