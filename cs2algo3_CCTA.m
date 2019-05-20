% Main Files to run examples
% This file computes the probability bound based on this new form of
% infinitesimal generator
clc; clear all; close all;
% Noise Level:
% gx:           0.5000    0.6000    0.7000  0.8000    0.9000    1.0000    1.1000    1.2000    1.3000    1.4000    1.5000
% alpha_st      0.9000    1.4000    1.5000  1.3000    1.3000    1.3000    1.3000    1.4000    1.5000    1.5000         0
% Monte Carlo:  0.2910    0.3964    0.4778  0.5520    0.6242    0.6710    0.7156    0.7628    0.7874    0.8322    0.8414
syms x1 x2
alpha = [.9];              % Range of alpha values for AB(x) =< -alpha*B(x) + Beta
gx = [.5];                  % Range of noise values
barrdeg =[10];
rho = 1;                    % Value over which we want to stay below
T = 2.0;                    % Time horizon (max time)
x0 = [-2 0];                % Initial condition
count = 1;                  % Counting variable
maxcount = 15;              % Max count value
cval = zeros(1, maxcount);  % Initialialize cmultiplier
beta = zeros(1, maxcount);  % Inialize beta storing
Pt = .10;

secondflag = 1;
for kk = 1:length(barrdeg)
    deg = barrdeg(kk);       % Degree of Polynomials we are searching
    
    % Pre-allocate optimate values
    probvalue = ones(length(alpha), length(gx));    % Initialize probability as ones
    alpha_st = zeros(length(alpha), length(gx));    % Initialize alpha as zero
    beta_st = zeros(length(alpha), length(gx));     % Initialize beta as zero
    gam_st = zeros(length(alpha),length(gx));       % Initialize gamma value
    Qtrace = [];
    probQ = [];
    beta = [];
    cvalstore = [];
    Qstore = [];
    for ii = 1:length(gx)
        for jj = 1:length(alpha)
            while count < maxcount
                if count == 1
                    u = 0;  % Control Input
                    [returnprob, retalpha, retbeta, retgam, Bx, prog] = cs2_verifctrl(alpha(jj), gx(ii), rho, T, x0, u, deg)
                    if probvalue(jj, ii) > returnprob && returnprob >= 0 && prog.solinfo.info.pinf == 0 && prog.solinfo.info.dinf == 0 
                       clc
                       probvalue(jj, ii) = returnprob
                       alpha_st(jj, ii) = retalpha
                       beta_st(jj, ii) = retbeta
                       gam_st(jj, ii) = retgam
                       Barriers(ii) = Bx;
                       beta(end+1) = retbeta;
                    elseif prog.solinfo.info.pinf == 1 || prog.solinfo.info.dinf == 1 
                        disp('Primal or dual infeasible problem.');
                    end
                else
                    if returnprob > Pt
                        beta(end+1) = .1*retbeta;
                        if secondflag
                            cfloor = 0;
                            secondflag = 0;
                        else
                            cfloor = 1.2*cval;
                        end
                    else 
                        beta(end+1) = 1.15*retbeta;
                        if secondflag
                            cfloor = 0;
                            secondflag = 0;
                        else
                            cfloor = .9*cval;
                        end
                    end
                    [u, traceQ, cval, Q] = cs2_initux(2, Bx, alpha_st(jj, ii), beta(end), gx(ii),cfloor)
                    [returnprob, retalpha, retbeta, retgam, Bx, prog] = cs2_verifctrl(alpha(jj), gx(ii), rho, T, x0, u, deg)
                     if returnprob >= 0 && prog.solinfo.info.pinf == 0 && prog.solinfo.info.dinf == 0 
                       probvalue(jj, ii) = returnprob
                       alpha_st(jj, ii) = retalpha
                       beta_st(jj, ii) = retbeta
                       gam_st(jj, ii) = retgam
                       Barriers(ii) = Bx;
                       if probvalue(jj, ii) < Pt
                            cvalstore(end+1) = cval;
                            Qtrace(end+1) = traceQ;
                            probQ(end+1) = returnprob;
                            Qstore(:,:,end+1) = Q;
                       end
                    elseif prog.solinfo.info.pinf == 1 || prog.solinfo.info.dinf == 1 
                        disp('Primal or dual infeasible problem.');
                    end

%                     if probvalue(jj, ii) < Pt
%                         break;
%                     end
                end
                    count = count + 1; 
            end
        end
    end
    if numel(probQ) > 0
        save(['quadctrl_polydeg' + string(deg) + '_' +  string(datetime('now','format','DDMMMyy_HHmmss'))])
    end
end
%%
figure;
plot(gx,probvalue, 'LineWidth', 2)
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
x1val = [-3:.05:3];
x2val = [-3:.05:3];

[X1 X2] = meshgrid(x1val, x2val);
x1 = X1;
x2 = X2;
Bvals = subs(Barriers(1));
figure;
[cm,c]=contour(X1, X2, Bvals, [0 1 3 20],'b-.','LineWidth',1.2);
clabel(cm, c,'Color','k','FontSize', 14, 'Interpreter','latex','LabelSpacing', 210);
grid on
grid minor
xlabel('$x_1$','Interpreter', 'latex')
ylabel('$x_2$','Interpreter', 'latex')
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')

patch([-4 4 4 -4],[2.25 2.25 4 4], 'red', 'FaceAlpha', .2,'EdgeColor','black')
lim = 3.0;
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
plot(0,0,'r*','Linewidth',2)
xlabel('$x_1$','Interpreter', 'latex')
ylabel('$x_2$','Interpreter', 'latex')
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')



% Plot Initial set
% xc1 = 0;
% xc2 = 0;
% r = .25;
% hold on;
% th = 0:pi/50:2*pi;
% x1unit = r * cos(th) + xc1;
% x2unit = r * sin(th) + xc2;
% plot(x1unit, x2unit,'g--', 'LineWidth', 2);
% xlabel('$x_1$','Interpreter', 'latex')
% ylabel('$x_2$','Interpreter', 'latex')
% set(gcf,'color','w');
% set(gca,'TickLabelInterpreter','latex')

% xlim([-5 5]);
% ylim([-5 5]);