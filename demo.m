% Minimal example for A*x - B*abs(x) = b.
% Run this script from any working directory. No files are generated.
packageRoot = fileparts(mfilename('fullpath'));
originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(fullfile(packageRoot, 'algos'));

rng(7, 'twister');
n = 32;
[Q, ~] = qr(randn(n));
A = Q * diag(linspace(2, 3, n)) * Q';
B = 0.2 * eye(n);
xTrue = randn(n, 1);
b = A * xTrue - B * abs(xTrue);

opts = struct('x0', zeros(n, 1), 'maxIter', 1000, ...
    'tolRes', 1e-10, 'storeHist', false, 'verbose', 0);
out = solve_gave_pgd_newton(A, B, b, opts);
relativeError = norm(out.x - xTrue) / max(1, norm(xTrue));

fprintf('GAVE-PGDN example (n = %d)\n', n);
fprintf('Iterations: %d\n', out.iter);
fprintf('Relative residual: %.3e\n', out.relRes);
fprintf('Relative solution error: %.3e\n', relativeError);
fprintf('Accepted Newton refinements: %d\n', out.nbNewtonAcc);
fprintf('Termination flag: %d\n', out.flag);
assert(out.relRes <= opts.tolRes, 'Demo failed to meet the residual tolerance.');
assert(relativeError <= 1e-8, 'Demo failed to recover the known solution.');
clear pathCleanup;
