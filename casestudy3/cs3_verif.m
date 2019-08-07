function [probvalue, alpha, betaval, Bxpolys, gam, prog] = cs3_verif(alpha, gx, rho, N, x0, u, deg)
    % Declare symbolic variables
    syms z x1 x2 betasym gamsym real
    % "z" is the noise variable for which we compute the moments
    EXP = 0;
    if ~exist('deg','var')
        % Check if degree if polynomial desired is given, if not default
        % to quartic
        deg = 4;
    end 
    betacoeff = .5;
    pJ = .95;
    pA = .5;

    % Setup Solver
    solver_opt.solver = 'sdpt3';
    fx = [  betacoeff*x1; ...
            pJ*x1 + pA*x2; ...
            ];            % Diffusion terms
    % Define program
    prog = sosprogram([x1, x2], [betasym, gamsym]);
    Zmon = monomials([x1, x2], [0:deg]);
    [prog, B] = sospolyvar(prog, Zmon, 'wscoeff');
    [prog, sig_u] = sospolyvar(prog, Zmon);
    [prog, sig_x] = sospolyvar(prog, Zmon);
    [prog, sig_o] = sospolyvar(prog, Zmon);
%     [prog, sig_1] = sospolyvar(prog, Zmon);
%     [prog, sig_2] = sospolyvar(prog, Zmon);
    % Define Inequalities
    prog = sosineq(prog, betasym);
%     prog = sosineq(prog, sig_1);
%     prog = sosineq(prog, sig_2);
    prog = sosineq(prog, sig_u);
    prog = sosineq(prog, sig_x);
    prog = sosineq(prog, sig_o);
    prog = sosineq(prog, B);
%     prog = sosineq(prog, B - sig_1*(x1) - sig_2*(x2));
    prog = sosineq(prog, gamsym);
    prog = sosineq(prog, 1 - gamsym - 1e-6);
    prog = sosineq(prog, B - sig_u*(x1^2 + x2^2 - 2^2 - 1e-3) - 1);
%     prog = sosineq(prog, B - sig_u*(x1 - 2 - 1e-3) - 1);
    prog = sosineq(prog, -B - sig_o*(1.5^2 - x1^2 - x2^2) + gamsym);
%     prog = sosineq(prog, -B - sig_o*(1^2 - (x1 - 1)^2 - (x2 - 1)^2) + gamsym);
    stdvar = gx;
    x1 = fx(1);
    x2 = fx(2) + z;

    Bsub = expand(subs(B));
    clear x1 x2;
    syms x1 x2 real;
    termlist = children(Bsub);

    % The following equations are for computer the Expected value
    for ii = 1:length(termlist)
        zcount = 0;
        x1count = 0;
        x2count = 0;

        % Initialize expectation as zero value
        EXPz = 0;
        
        factoredterm = factor(termlist(ii));
        % Count power of each variable in monomial
        for jj = 1:length(factoredterm)
            
            % Count the power of the "noise" term
            if isequaln(factoredterm(jj),z)
                zcount = zcount + 1;
            end
            
            % Count nubmer of states
            if isequaln(factoredterm(jj),x1)
                x1count = x1count + 1;
            end
            
            if isequaln(factoredterm(jj),x2)
                x2count = x2count + 1;
            end
        end
        
        % Check to see if symbolic term has no noise variables
        % May or may not have state variables
        if zcount == 0
            EXPz = EXPz + termlist(ii);
        end
        
        % Check if there are "noise" terms but not x variables
        if zcount > 0 && x1count == 0 && x2count == 0
            if mod(zcount,2) == 1
                EXPz = 0;
            elseif mod(zcount,2) == 0
                EXPz = prod(factoredterm(find(factoredterm~=z)))*prod([1:2:zcount])*stdvar^zcount;
            end
        end
        
        % Check if there are "noise" terms and x variables
        if zcount > 0 && (x1count > 0 || x2count > 0)
            if mod(zcount,2) == 1
                EXPz = 0;
            elseif mod(zcount,2) == 0
                EXPz = prod(factoredterm(find(factoredterm~=z)))*prod([1:2:zcount])*stdvar^zcount;
            end
        end
        
        EXP = EXP + EXPz;
    end
    
%   Expectation evolution
    prog = sosineq(prog, -EXP + B/alpha + betasym - sig_x*(2^2 - x1^2 - x2^2));
%     prog = sosineq(prog, -EXP + B/alpha + betasym - sig_x*(2 - x1));

    if alpha == 1
        objfunc = gamsym + betasym;
        prog = sossetobj(prog, objfunc);
    else
        objfunc = gamsym + betasym;
%         objfunc =  betasym;
        prog = sossetobj(prog, objfunc);
    end

    % Solve Program 
    prog = sossolve(prog, solver_opt);
    Bxpolys = sosgetsol(prog, B);
    betaval = double(sosgetsol(prog, betasym));
    gam = sosgetsol(prog, gamsym);
    x1 = x0(1);
    x2 = x0(2);

    if alpha == 1
        probvalue = (subs(Bxpolys) + betaval*N);
    elseif betaval*alpha/(alpha - 1) <= 1 && alpha > 1
        picomp = 1;
        for ii = 1:N-1
            picomp = picomp*(1 - betaval);
        end
        probvalue = 1 - ( 1 - subs(Bxpolys))*picomp;
    elseif betaval*alpha/(alpha - 1) > 1 && alpha > 1
        probvalue = subs(Bxpolys)*alpha^(-N) + (1 - alpha^-N)*alpha*betaval/(alpha - 1);
    end
    
    probvalue = double(probvalue);
end

