function [probvalue, alpha, betaval, gamval, Bxpolys, prog] = cs2_verifctrl(alpha, gx, rho, T, x0, u, deg)
    % Declare symbolic variables
    syms x1 x2 betasym gamsym real
    if ~exist('deg','var')
        % Check if degree if polynomial desired is given, if not default
        % to quartic
        deg = 4;
    end 

    % Setup Solver
    solver_opt.solver = 'sdpt3';
    fx = [  x2; ...
        -x1 - x2 - .5*x1^3] + [0; 1]*u;   
    
    % Define program
    prog = sosprogram([x1 x2], [betasym, gamsym]);
    Zmon = monomials([x1 x2], [0:deg]);
    
    % Define Polynomials
    [prog, B] = sospolyvar(prog, Zmon, 'wscoeff');  % Barrier
    [prog, sig_u] = sospolyvar(prog, Zmon);         % SOS unsafe poly
    [prog, sig_s] = sospolyvar(prog, Zmon);         % SOS initial state poly
    [prog, sig_o] = sospolyvar(prog, Zmon);         % SOS initial state poly

    % Define generator
    dBdx1 = diff(B, x1);
    dBdx2 = diff(B, x2);
    d2Bd2x11 = diff(dBdx1, x1);
    d2Bd2x22 = diff(dBdx2, x2);
    d2Bd2x12 = diff(dBdx1, x2);
    d2Bd2x21 = diff(dBdx2, x1);
    dBdx = [d2Bd2x11, d2Bd2x12; d2Bd2x21, d2Bd2x22];
    AB = [dBdx1, dBdx2]*fx + .5*trace(gx'*dBdx*gx);

    % Define Inequalities
    prog = sosineq(prog, sig_u);    % Unsafe set multiplier
    prog = sosineq(prog, sig_s);    % Safe set multiplier
    prog = sosineq(prog, sig_o);    % Inital set constraint
    prog = sosineq(prog, betasym);
%     prog = sosineq(prog, prog.symdecvartable(2));
%     prog = sosineq(prog, 1 - prog.symdecvartable(2) - 1e-6);
    prog = sosineq(prog, gamsym);
    prog = sosineq(prog, 1 - gamsym - 1e-6);
%     prog = sosineq(prog, alpha - betasym*T);
   
    % Encode constraints on sets
    prog = sosineq(prog, B);        % Declare barrier as non-negative
    prog = sosineq(prog, B - sig_u*(x2 - 2.25) - 1);
    prog = sosineq(prog, -alpha*B + betasym - AB - sig_s*(2.25 - x2));
    prog = sosineq(prog, -B - sig_o*(.1^2 - (x1 + 2)^2 - x2^2) + gamsym); % Initial set constraint ()
%     prog = sosineq(prog, betasym - alpha*B) ; % Initial set constraint ()
    
    if alpha == 0
        objfunc = prog.symdecvartable(2) + betasym;
        prog = sossetobj(prog, objfunc);        
    else
%         objfunc = betasym*T;
        objfunc = prog.symdecvartable(2) + betasym;
        prog = sossetobj(prog, objfunc);
    end

    % Solve Program 
    prog = sossolve(prog, solver_opt);
    Bxpolys = sosgetsol(prog, B);
    betaval = double(sosgetsol(prog, betasym));
    gamval = double(sosgetsol(prog, gamsym));
%     gamval = 0;
    x1 = x0(1);
    x2 = x0(2);
    
    if alpha == 0
        probvalue = (subs(Bxpolys) + betaval*T)/rho;
    elseif rho >= betaval/alpha
        probvalue = 1 - ( 1 - subs(Bxpolys)/rho)*exp(-betaval*T) ;
    elseif rho < betaval/alpha
        probvalue =  (subs(Bxpolys) + (exp(betaval*T) - 1)*(betaval/alpha) )/(rho*exp(betaval*T)) ;
    end
    
    probvalue = double(probvalue);
end

