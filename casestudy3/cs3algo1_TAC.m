%% Probability Bound
clear all; clc; close all;

alpha_ = [1:5:20];   % Range of alpha values for AB(x) =< -alpha*B(x) + Beta
gx = [0:.1:.3];              % Range of noise values

rho = 1;                % Value over which we want to stay below
N = 2;                % Time horizon (max time)
x0 = [1; 1];                 % Initial condition
deg = 8;                % Degree of Polynomials we are searching

% Pre-allocate optimate values
probvalue = ones(1,length(gx));
gamprob = ones(1,length(gx));
alpha_st = zeros(1,length(gx));
beta_st = zeros(1,length(gx));
gam_st = zeros(1,length(gx));
u_st = zeros(1,length(gx));

for ii = 1:length(gx)
    for jj = 1:length(alpha_)
        [returnprob, retalpha, retbeta, Bx, retgam, prog] = cs3_verif(alpha_(jj), gx(ii), rho, N, x0, 0, deg)
        if probvalue(ii) > returnprob && returnprob >= 0 && prog.solinfo.info.pinf == 0 && prog.solinfo.info.dinf == 0
           clc;
           probvalue(ii) = returnprob
           alpha_st(ii) = retalpha
           beta_st(ii) = retbeta
           gam_st(ii) = retgam
           Barriers(ii) = Bx;
        end
    end
end

% save(["probboundrefined_deg"+ string(deg) + "_" + datestr(now,'mmddyy_HHMMSS') + ".mat"])
%%
figure;
plot(gx, probvalue,'LineWidth',2,'DisplayName',["Upper bound"]);
grid on
grid minor
xlabel('$\sigma(x)$','Interpreter','latex');
ylabel('Probability','Interpreter','latex');
l = legend('show');
set(l,'Interpreter','latex');
set(gca,'FontSize',14,'TickLabelInterpreter','latex');
set(gcf,'color','w');

savefig('probbound');
lim = 4.0;
x1 = [-lim:.5:lim];
x2 = [-lim:.5:lim];

[X1, X2] = meshgrid(x1,x2);
x1 = X1;
x2 = X2;
Bvals = double(subs(Barriers(1)));
figure;
[c, h] = contour(X1,X2,Bvals,'b-.','DisplayName','$B(x)$ Level Sets','LineWidth',1.2);
hold on;
clabel(c, h,'Color','k','FontSize', 14, 'Interpreter','latex','LabelSpacing', 210);

x1 = 0;
x2 = 0;
r = 2;
hold on;
xlim([-lim lim])
ylim([-lim lim])
th = 0:pi/50:2*pi;
x1unit = r * cos(th) + x1;
x2unit = r * sin(th) + x2;
plot(x1unit, x2unit,'g--','LineWidth',2,'DisplayName','$\mathcal{X}_0$');
patch(x1unit, x2unit, 'white', 'FaceAlpha', 0,'EdgeColor','none','HandleVisibility','off')
patch([-lim r x1unit r -lim],[-lim -lim x2unit lim lim], 'red', 'FaceAlpha', .2,'EdgeColor','none','HandleVisibility','off')
patch([r lim lim r],[-lim -lim lim lim], 'red', 'FaceAlpha', .2,'EdgeColor','none','HandleVisibility','off')
% patch([2 lim lim 2],[0 0 lim lim], 'red', 'FaceAlpha', .2,'EdgeColor','none','HandleVisibility','off')
text(2.5,3.5,'$\mathcal{X}_u$','Interpreter','latex','FontSize',20);

grid on
grid minor
xlabel('$x_1$','Interpreter', 'latex')
ylabel('$x_2$','Interpreter', 'latex')
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex','FontSize', 14)
