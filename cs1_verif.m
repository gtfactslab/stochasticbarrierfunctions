function [probvalue, alpha, betaval, Bxpolys, gam, prog] = cs1_verif(alpha, gx, rho, T, x0, u, deg)
    % Declare symbolic variables
    syms x betasym gamsym real

    if ~exist('deg','var')
        % Check if degree if polynomial desired is given, if not default
        % to quartic
        deg = 4;
    end 

    % Setup Solver
    solver_opt.solver = 'sdpt3';
    fx = -1*x + u;            % Diffusion terms

    % Define program
    prog = sosprogram([x], [betasym, gamsym]);
    Zmon = monomials([x], [0:deg]);
    [prog, B] = sospolyvar(prog, Zmon, 'wscoeff');
    [prog, sig_u] = sospolyvar(prog, Zmon);
    [prog, sig_x] = sospolyvar(prog, Zmon);
    [prog, sig_o] = sospolyvar(prog, Zmon);

    % Define generator
    dBdx = diff(B, x);
    d2Bdx2 = diff(dBdx, x);
    AB = diff(B, x)*fx + .5*trace(gx'*d2Bdx2*gx);

    % Define Constraints
    prog = sosineq(prog, betasym);
    prog = sosineq(prog, sig_u);
    prog = sosineq(prog, sig_x);
    prog = sosineq(prog, sig_o);
    prog = sosineq(prog, B);
    prog = sosineq(prog, gamsym);
    prog = sosineq(prog, 1 - gamsym - 1e-6);
    prog = sosineq(prog, B - sig_u*(x.^2 - 1) - 1);
    prog = sosineq(prog, -B - sig_o*(.2^2 - x.^2) + gamsym);
    prog = sosineq(prog, -AB -alpha*B + betasym - sig_x*(1 - 1e-3 - x.^2));
    
    if alpha == 0
        objfunc = prog.symdecvartable(2) + betasym;
        prog = sossetobj(prog, objfunc);

    else
        objfunc = gamsym + betasym;
        prog = sossetobj(prog, objfunc);
    end

    % Solve Program 
    prog = sossolve(prog, solver_opt);
    Bxpolys = sosgetsol(prog, B);
    betaval = double(sosgetsol(prog, betasym));
%       gam = 1;
    gam = sosgetsol(prog, gamsym);
    x= x0;
    if alpha == 0
        probvalue = (subs(Bxpolys) + betaval*T)/rho;
    elseif rho >= betaval/alpha
        probvalue = 1 - ( 1 - subs(Bxpolys)/rho)*exp(-betaval*T) ;
    elseif rho < betaval/alpha
        probvalue =  (subs(Bxpolys) + (exp(betaval*T) - 1)*(betaval/alpha) )/(rho*exp(betaval*T)) ;
    end
    
    probvalue = double(probvalue);
end

