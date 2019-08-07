%% Initialize u(x)
function [ux, traceQ, cval, Qvals] = cs4_initux(deg, B, alpha, beta, gx, cfloor)
    syms x1 x2 c z real;
    solver_opt.solver = 'sdpt3';
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
    prog = sosprogram([x1 x2], [c]);
    Zmon = monomials([x1 x2],[0:deg]);

    dim = length(Zmon);
    [prog, Q] = sospolymatrixvar(prog, monomials([x1],0),[dim, dim]); 

    ux = Zmon'*Q*Zmon;
    onevec = ones(numel(Q), numel(Q));
    [prog, sig_x] = sospolyvar(prog, Zmon);
    
    stdvar = gx;
    x1 = fx(1)+ ux;
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
                EXPz = prod([1:2:zcount])*stdvar^zcount;
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
    
    prog = sosineq(prog, -EXP + B/alpha + beta - sig_x*(2 - x1));
    
    prog = sosineq(prog, c - cfloor);
    for ii = 1:numel(Q)
        prog = sosineq(prog, c - Q(ii));
        prog = sosineq(prog, Q(ii) + c);
    end
    objfunc = c;
    prog = sossetobj(prog, objfunc);

    prog = sossolve(prog, solver_opt);

    Qvals = double(sosgetsol(prog, Q));
    traceQ = trace(Qvals);
    ux = Zmon'*Qvals*Zmon;
    cval = sosgetsol(prog, c);
    
end

