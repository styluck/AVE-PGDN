# GAVE-PGDN

Minimal MATLAB implementation of projected gradient descent with adaptive Newton refinement for the generalized absolute value equation

```text
A*x - B*abs(x) = b
```

This package contains the solver copied unchanged from the current research code and one small synthetic example. It does not contain the paper's benchmark experiments, datasets, plots, or saved results. The demo illustrates usage rather than reproducing the numerical study in the paper.

## Quick start

Open MATLAB in this directory and run:

```matlab
run('demo.m')
```

Alternatively, run the demo by its absolute path. The demo sets up and restores the MATLAB search path automatically, uses a fixed random seed, and prints iteration count, relative residual, solution error, and accepted Newton refinements. It creates no output files. The matrix A is positive definite with minimum eigenvalue 2 and B = 0.2*I; the right-hand side is constructed from a known solution.

## Requirements

MATLAB; no additional toolboxes or third-party code are required. Validation details are in VALIDATION.md. The rectangular Newton branch uses the built-in function `lsqminnorm`; use a MATLAB version providing that function when solving rectangular problems.

## Solver usage

```matlab
opts = struct('x0', zeros(size(A,2),1), 'maxIter', 1000, 'tolRes', 1e-10);
out = solve_gave_pgd_newton(A, B, b, opts);
x = out.x;
```

For direct calls, first add `algos` to the MATLAB path. A and B must have identical dimensions, and b must have size(A,1) entries.

| Option | Default | Meaning |
| --- | --- | --- |
| x0 | zero vector | Initial solution estimate |
| maxIter | 5000 | Maximum projected-gradient iterations |
| tolRes | 1e-10 | Relative residual tolerance |
| tolStep | 1e-12 | Projected step stopping tolerance |
| stepType | 'bb' | 'bb' or 'fixed' step selection |
| useLineSearch | true | Enable backtracking |
| useNewton | true | Enable adaptive linearized refinement |
| newtonTol | 1e-3 | Residual threshold for attempting refinement |
| signStableIts | 5 | Sign stability threshold for attempting refinement |
| storeHist | false | Return iteration history in out.hist |
| verbose | 0 | Print iteration diagnostics |

See the solver source for all options. BB steps and adaptive refinement are implementation choices; this demo does not establish the assumptions of any particular convergence theorem.

Important outputs are `out.x`, `out.iter`, `out.relRes`, `out.compInf`, `out.nbNewtonTry`, and `out.nbNewtonAcc`. Here `out.relRes = norm(A*x-B*abs(x)-b)/max(1,norm(b))`. `out.nbLP` counts projected-gradient iterations, not linear-programming solves.

Termination flags: 0 = iteration limit; 1 = residual tolerance met; 2 = small projected step without meeting residual tolerance; 3 = line search failed. Always check the final residual as well as the flag.

## Contents

- `algos/solve_gave_pgd_newton.m`: self-contained GAVE solver, including local helper functions.
- `demo.m`: fixed-seed, small square GAVE example with a known solution.
- `VALIDATION.md`: checks performed on this package.

## Paper and license

Before public release, add the final paper title, author list, and citation/link here. A distribution license has not yet been selected; no license grant is supplied by this package. Add the copyright holder's approved LICENSE before publishing as licensed open-source software.
