# Validation

Validated on 2026-10-05 using MATLAB R2024b on Windows.

- Ran demo.m from the project workspace and by absolute path from MATLAB's temporary directory. Both runs completed successfully.
- Confirmed that the demo restores the MATLAB search path.
- The demo asserts relative residual <= 1e-10 and relative solution error <= 1e-8.
- Observed: 12 iterations, relative residual 2.759e-16, relative solution error 2.724e-16, one accepted Newton refinement, termination flag 1.
- Verified the released solver is byte-for-byte identical to _paper/code/algos/solve_gave_pgd_newton.m using SHA-256: 9BCE567D38B605AEC6099A644AC46F1D208DA8FD7659DEC231CE34DB416BEC67.

This validates the small square example only; the rectangular branch and the paper's full experiments were not run for this package.
