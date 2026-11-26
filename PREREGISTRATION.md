# Plan written before the results

Date: 2026-11-23

The model is the transverse field Ising chain, H = -sum X_i X_{i+1} - h sum Z_i, with J = 1. The central charge of its critical point at h = 1 is known exactly: c = 1/2. This project measures c from the ground state entanglement entropy and compares it with 1/2. It reproduces a known result. It does not discover one.

I fixed these choices before running anything:

1. Main estimate: c from the half chain entropy of the open chain at h = 1, solved with free fermions at L = 1024 and L = 2048, using c = 6 (S(2048) - S(1024)) / ln 2.
2. Its uncertainty: the change in that estimate from the previous doubling, |c(1024->2048) - c(512->1024)|. There is no statistical noise in this project, because every number is exact up to rounding. The only error is from the finite chain length.
3. Test: I will call the result consistent with 1/2 if the difference is less than 3 times that uncertainty. If it is not, I will report that and will not change the estimator.
4. Supporting estimates: a fit of S(l) against the Calabrese and Cardy form for the open chain at L = 1024 (l = 64 to 512, step 16), a fit for the periodic chain at L = 20 by exact diagonalization (l = 2 to 18), and c from the periodic half chain entropy at L = 18 and 20.
5. Checks that must pass first: the exact diagonalization matches the exact periodic energy formula, and the free fermion solution matches exact diagonalization for open chains up to L = 20.
