%% Initialize u(x)
function [ux, traceQ, cval, Qvals] = cs2_initux(deg, B, alpha, beta, gx, cfloor)
    syms x1 x2 c real;
    mu = 1;
    solver_opt.solver = 'sdpt3';

    prog = sosprogram([x1 x2], [c]);
    Zmon = monomials([x1 x2],[0:deg]);

    dim = length(Zmon);
    [prog, Q] = sospolymatrixvar(prog, monomials([x1],0),[dim, dim]); 

    ux = Zmon'*Q*Zmon;
    onevec = ones(numel(Q), numel(Q));
    [prog, sig_s] = sospolyvar(prog, Zmon);         % SOS initial state poly

    fx = [  x2; ...
        -x1 - x2 - .5*x1^3]  + [0; 1]*ux;

    dBdx1 = diff(B, x1);
    dBdx2 = diff(B, x2);
    d2Bd2x11 = diff(dBdx1, x1);
    d2Bd2x22 = diff(dBdx2, x2);
    d2Bd2x12 = diff(dBdx1, x2);
    d2Bd2x21 = diff(dBdx2, x1);
    dBdx = [d2Bd2x11, d2Bd2x12; d2Bd2x21, d2Bd2x22];
    AB = [dBdx1, dBdx2]*fx + .5*trace(gx'*dBdx*gx);
    
    prog = sosineq(prog, -alpha*B + beta - AB - sig_s*(2.25 - x2));
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

