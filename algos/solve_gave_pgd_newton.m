function out = solve_gave_pgd_newton(A, B, b, opts)
%SOLVE_GAVE_PGD_NEWTON Hybrid BB-PGD + adaptive linearized refinement for GAVE.

    if nargin < 4
        opts = struct();
    end

    [m, n] = size(A);
    if ~isequal(size(B), [m, n])
        error('A and B must have the same size.');
    end
    b = b(:);
    if numel(b) ~= m
        error('Dimension mismatch: length(b) must equal size(A,1).');
    end

    maxIter       = get_opt(opts, 'maxIter', 5000);
    tolRes        = get_opt(opts, 'tolRes', 1e-10);
    tolStep       = get_opt(opts, 'tolStep', 1e-12);
    verbose       = get_opt(opts, 'verbose', 0);
    storeHist     = get_opt(opts, 'storeHist', false);

    alpha0        = get_opt(opts, 'alpha0', 1.0);
    alphaMin      = get_opt(opts, 'alphaMin', 1e-12);
    alphaMax      = get_opt(opts, 'alphaMax', 1e12);
    beta          = get_opt(opts, 'beta', 0.5);
    sigma         = get_opt(opts, 'sigma', 1e-4);
    useLineSearch = get_opt(opts, 'useLineSearch', true);
    stepType      = lower(get_opt(opts, 'stepType', 'bb'));
    bbType        = upper(get_opt(opts, 'bbType', 'ALT'));

    useNewton     = get_opt(opts, 'useNewton', true);
    newtonTol     = get_opt(opts, 'newtonTol', 1e-3);
    signStableIts = get_opt(opts, 'signStableIts', 5);
    newtonEta     = get_opt(opts, 'newtonEta', 0.8);
    maxNewton     = get_opt(opts, 'maxNewton', 3);
    newtonDamped  = get_opt(opts, 'newtonDamped', true);
    newtonGammas  = get_opt(opts, 'newtonGammas', [1, 0.5, 0.25, 0.1]);

    Mplus = A - B;
    Mminus = -A - B;

    has_xp0 = isfield(opts, 'xp0');
    has_xm0 = isfield(opts, 'xm0');
    has_x0 = isfield(opts, 'x0');

    if has_xp0 && has_xm0
        xp = max(opts.xp0(:), 0);
        xm = max(opts.xm0(:), 0);
        [xp, xm] = proj_comp_pair(xp, xm);
    elseif has_x0
        x0 = opts.x0(:);
        if numel(x0) ~= n
            error('opts.x0 has incompatible dimension.');
        end
        xp = max(x0, 0);
        xm = max(-x0, 0);
    else
        xp = zeros(n, 1);
        xm = zeros(n, 1);
    end

    if storeHist
        hist.iter = [];
        hist.f = [];
        hist.res = [];
        hist.relRes = [];
        hist.step = [];
        hist.alpha = [];
        hist.gradInf = [];
        hist.signChanges = [];
        hist.newtonTried = [];
        hist.newtonTaken = [];
    end

    r = Mplus * xp + Mminus * xm - b;
    f = 0.5 * (r' * r);
    gxp = Mplus' * r;
    gxm = Mminus' * r;

    z_prev = [];
    g_prev = [];
    alpha_prev = alpha0;

    sgn_prev = sign_pattern(xp - xm);
    stableCount = 0;

    nbLinSys = 0;
    nbLP = 0;
    nbNewtonTry = 0;
    nbNewtonAcc = 0;
    flag = 0;

    for k = 1:maxIter
        x = xp - xm;
        Fx = A * x - B * abs(x) - b;
        relRes = norm(Fx) / max(1, norm(b));

        if relRes <= tolRes
            flag = 1;
            nbLP = k - 1;
            break;
        end

        if strcmp(stepType, 'bb')
            if k == 1 || isempty(z_prev)
                alpha = alpha0;
            else
                z = [xp; xm];
                g = [gxp; gxm];
                s = z - z_prev;
                y = g - g_prev;

                sty = s' * y;
                yty = y' * y;
                sts = s' * s;

                switch bbType
                    case 'BB1'
                        if sty > 0
                            alpha = sts / sty;
                        else
                            alpha = alpha_prev;
                        end
                    case 'BB2'
                        if yty > 0
                            alpha = sty / yty;
                        else
                            alpha = alpha_prev;
                        end
                    otherwise
                        if mod(k, 2) == 0
                            if sty > 0
                                alpha = sts / sty;
                            else
                                alpha = alpha_prev;
                            end
                        else
                            if yty > 0
                                alpha = sty / yty;
                            else
                                alpha = alpha_prev;
                            end
                        end
                end

                if ~isfinite(alpha) || alpha <= 0
                    alpha = alpha_prev;
                end
            end
        elseif strcmp(stepType, 'fixed')
            alpha = alpha0;
        else
            error('Unknown opts.stepType. Use ''bb'' or ''fixed''.');
        end

        alpha = min(max(alpha, alphaMin), alphaMax);

        z_curr = [xp; xm];
        g_curr = [gxp; gxm];

        accepted = false;
        alpha_trial = alpha;

        while true
            yxp = xp - alpha_trial * gxp;
            yxm = xm - alpha_trial * gxm;

            [xp_bar, xm_bar] = proj_comp_pair(yxp, yxm);
            x_bar = xp_bar - xm_bar;

            r_bar = Mplus * xp_bar + Mminus * xm_bar - b;
            f_bar = 0.5 * (r_bar' * r_bar);
            stepNorm = norm([xp_bar - xp; xm_bar - xm], inf);

            if ~useLineSearch
                accepted = true;
                break;
            end

            dxp = xp_bar - xp;
            dxm = xm_bar - xm;
            rhs = f - (sigma / alpha_trial) * (norm(dxp)^2 + norm(dxm)^2);

            if f_bar <= rhs
                accepted = true;
                break;
            end

            alpha_trial = beta * alpha_trial;
            if alpha_trial < alphaMin
                break;
            end
        end

        if ~accepted
            flag = 3;
            nbLP = k - 1;
            break;
        end

        F_bar = A * x_bar - B * abs(x_bar) - b;
        relRes_bar = norm(F_bar) / max(1, norm(b));

        tookNewton = false;
        triedNewton = false;

        sgn_bar = sign_pattern(x_bar);
        if isequal(sgn_bar, sgn_prev)
            stableCount = stableCount + 1;
        else
            stableCount = 0;
        end

        triggerNewton = useNewton && (relRes_bar <= newtonTol || stableCount >= signStableIts);

        if triggerNewton
            triedNewton = true;
            nbNewtonTry = nbNewtonTry + 1;

            x_best = x_bar;
            rel_best = relRes_bar;
            xN = x_bar;

            for j = 1:maxNewton
                zN = linearized_solve(A, B, b, xN);
                if isempty(zN) || any(~isfinite(zN))
                    break;
                end

                nbLinSys = nbLinSys + 1;

                if newtonDamped
                    acceptedNewtonInner = false;
                    for gamma = newtonGammas
                        xCand = xN + gamma * (zN - xN);
                        relCand = norm(A * xCand - B * abs(xCand) - b) / max(1, norm(b));
                        if relCand <= newtonEta * rel_best
                            xN = xCand;
                            rel_best = relCand;
                            x_best = xCand;
                            acceptedNewtonInner = true;
                            break;
                        end
                    end
                    if ~acceptedNewtonInner
                        break;
                    end
                else
                    relCand = norm(A * zN - B * abs(zN) - b) / max(1, norm(b));
                    if relCand <= newtonEta * rel_best  % relRes_bar %
                        xN = zN;
                        rel_best = relCand;
                        x_best = zN;
                    else
                        break;
                    end
                end

                if rel_best <= tolRes
                    break;
                end
            end

            if rel_best < relRes_bar
                x_new = x_best;
                xp_new = max(x_new, 0);
                xm_new = max(-x_new, 0);
                tookNewton = true;
                nbNewtonAcc = nbNewtonAcc + 1;
            else
                x_new = x_bar;
                xp_new = xp_bar;
                xm_new = xm_bar;
            end
        else
            x_new = x_bar;
            xp_new = xp_bar;
            xm_new = xm_bar;
        end

        xp = xp_new;
        xm = xm_new;
        x = x_new;

        r = Mplus * xp + Mminus * xm - b;
        f = 0.5 * (r' * r);
        gxp = Mplus' * r;
        gxm = Mminus' * r;

        z_prev = z_curr;
        g_prev = g_curr;
        alpha_prev = alpha_trial;
        sgn_prev = sign_pattern(x);

        nbLP = k;

        if storeHist
            hist.iter(end+1,1) = k;
            hist.f(end+1,1) = f;
            hist.res(end+1,1) = norm(r);
            hist.relRes(end+1,1) = norm(A * x - B * abs(x) - b) / max(1, norm(b));
            hist.step(end+1,1) = stepNorm;
            hist.alpha(end+1,1) = alpha_trial;
            hist.gradInf(end+1,1) = norm([gxp; gxm], inf);
            hist.signChanges(end+1,1) = sum(sign_pattern(x) ~= sgn_bar);
            hist.newtonTried(end+1,1) = triedNewton;
            hist.newtonTaken(end+1,1) = tookNewton;
        end

        if verbose && (k == 1 || mod(k, 50) == 0 || tookNewton || k == maxIter)
            fprintf(['GAVE-PGDN iter %5d: relRes = %.3e, step = %.3e, alpha = %.3e, ' ...
                     'newtonTried = %d, newtonTaken = %d\n'], ...
                     k, norm(A * x - B * abs(x) - b) / max(1, norm(b)), ...
                     stepNorm, alpha_trial, triedNewton, tookNewton);
        end

        if stepNorm <= tolStep && ~tookNewton
            if norm(A * x - B * abs(x) - b) / max(1, norm(b)) <= tolRes
                flag = 1;
            else
                flag = 2;
            end
            break;
        end
    end

    x = xp - xm;

    out = struct();
    out.x = x(:);
    out.xp = xp(:);
    out.xm = xm(:);
    out.iter = nbLP;
    out.flag = flag;
    out.relRes = norm(A * out.x - B * abs(out.x) - b) / max(1, norm(b));
    out.compInf = norm(out.xp .* out.xm, inf);
    out.nbLP = nbLP;
    out.nbLin = nbLinSys;
    out.nbNewtonTry = nbNewtonTry;
    out.nbNewtonAcc = nbNewtonAcc;

    if storeHist
        out.hist = hist;
    end
end

function z = linearized_solve(A, B, b, x)
    d = sign_diag_vector(x);
    J = A - B * spdiags(d, 0, numel(d), numel(d));
    z = [];

    try
        if size(J, 1) == size(J, 2)
            % Never densify a sparse Newton matrix merely to estimate its
            % reciprocal condition number.  The old rcond(full(J)) check
            % made the nominally sparse large-scale experiment require
            % O(n^2) memory before the sparse direct solve was attempted.
            if ~issparse(J)
                if rcond(J) < 1e-15
                    return;
                end
            end
            zcand = J \ b;
        else
            zcand = lsqminnorm(J, b);
        end
        if all(isfinite(zcand))
            z = zcand;
        end
    catch
        z = [];
    end
end

function [xp, xm] = proj_comp_pair(a, b)
    ap = max(a, 0);
    bp = max(b, 0);

    xp = zeros(size(ap));
    xm = zeros(size(bp));

    idx_u = (ap >= bp);
    idx_v = ~idx_u;

    xp(idx_u) = ap(idx_u);
    xm(idx_u) = 0;
    xp(idx_v) = 0;
    xm(idx_v) = bp(idx_v);
end

function s = sign_pattern(x)
    s = zeros(size(x));
    s(x > 0) = 1;
    s(x < 0) = -1;
end

function d = sign_diag_vector(x)
    d = -ones(size(x));
    d(x >= 0) = 1;
end

function val = get_opt(opts, name, defaultVal)
    if isfield(opts, name) && ~isempty(opts.(name))
        val = opts.(name);
    else
        val = defaultVal;
    end
end
