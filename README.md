# Random matrix theory-based method for evaluating vulnerability of physical unclonable functions to modeling attacks

This MATLAB code accompanies the article of the same title. Given paired challenges and responses from a physical unclonable function (PUF), it characterizes how response variance is distributed across modes, how much of the spectrum exceeds a random-matrix reference, how neighboring eigenvalues are spaced, and how strongly challenges and responses depend on each other. These diagnostics describe the measured dataset. Quantifying prediction accuracy requires a separate evaluation on held-out challenge–response pairs.

## Run the analysis

1. Place `main.m` and `puf_metrics.m` in the same folder and open that folder in MATLAB. If the downloaded filenames end in `(1).m`, rename them to match these names.
2. Prepare a challenge matrix `C` and a response matrix `R`. Each row must describe the same observation in both files. Both matrices must contain finite numeric or logical values and have at least two rows. Real and complex measurements are supported.
3. Set the input paths and pairing mode at the top of `main.m`. The default configuration is shown below.

```matlab
challengeFile = 'challenges.mat';    % Must contain a variable named C
responseFile  = 'responses.mat';     % Must contain a variable named R
pairing = 'rows';
```

Numeric CSV files are also supported. Set the paths to the corresponding CSV filenames and store only the matrix values, with observations in rows.

4. Run `main` in the MATLAB Command Window. The script displays a one-row table and leaves it in the workspace as `results`. To export the table, run `writetable(results,'metrics.csv')`.

## Choose the pairing mode

| Mode | Input dimensions | One observation for distance correlation |
| --- | --- | --- |
| `'rows'` | `C` is N × D and `R` is N × P. D and P may differ. | A complete challenge vector paired with its response vector. For example, seven control settings can be paired with ten measured output powers. |
| `'entries'` | Both matrices are N × P and share the same grid. | One scalar input paired with the scalar response at the same position. On a current–wavelength grid, if each column of `C` repeats the current vector, dCor measures current–response dependence pooled across wavelengths. |

The choice affects only `dCor`. The three spectral metrics always use the full response matrix. The supplied `main.m` explicitly selects `'rows'`. A direct call to `puf_metrics(C,R)` defaults to `'entries'`. Pass the mode explicitly when switching between datasets.

## Output metrics

For a response matrix with N rows and P columns, the function centers each column and scales each varying column to unit sample standard deviation. Constant columns are retained as zeros. If the standardized matrix has singular values σᵢ, the covariance eigenvalues are λᵢ = σᵢ²/(N − 1).

| Output column | Calculation and interpretation |
| --- | --- |
| `rho_eff` | Normalized effective rank. With p = min(N − 1, P) and weights wᵢ = λᵢ/Σλⱼ, the function computes exp(−Σ wᵢ ln wᵢ)/p over positive weights. A single mode carrying all the variance gives 1/p. Equal variance across all p modes gives 1. |
| `S_spike` | Fraction of total spectral energy in excess of the Marchenko–Pastur upper edge. The edge is λ₊ = (Σλᵢ/P)(1 + √(P/(N − 1)))². Each eigenvalue above the edge contributes only its excess λᵢ − λ₊. A value of zero means all eigenvalues lie at or below the edge. |
| `beta_Brody` | Maximum-likelihood Brody fit to unfolded nearest-neighbor spacings. Values near 0 indicate Poisson-like spacings. Values near 1 indicate level repulsion close to the Gaussian orthogonal ensemble reference. |
| `dCor` | Sample distance correlation between paired challenges and responses, computed on the original data using Euclidean distances. Higher values indicate stronger measured dependence, including nonlinear dependence. The estimator is biased and can give positive values for independent samples. |

For the Brody fit, the code retains singular values above 10⁻¹⁰ times the largest singular value, applies Gaussian unfolding to the corresponding covariance levels, and excludes the outer 10% of levels on each side (rounded down) when forming spacings. For dCor, the code uses all paired observations up to 1,000. Above that limit, it selects 1,000 paired observations without replacement using a reproducible seed determined by the dimensions of `R`.

## Quick check in MATLAB

This example runs directly on in-memory matrices. All response columns are proportional to the same input, so the standardized response has one mode.

```matlab
C = (1:101)';
R = C * (1:50);
T = puf_metrics(C,R,'rows');
disp(T)
```

The expected values are `rho_eff = 0.02` (1/50), `S_spike ≈ 0.9417`, and `dCor = 1` up to rounding. `beta_Brody` is `NaN` because one eigenvalue supplies no adjacent spacings.

## When a result is NaN

- All three spectral metrics are `NaN` if the centered response has zero total energy, as happens when every response feature is constant across observations.
- `beta_Brody` is also `NaN` if fewer than 30 interior spacings remain, the unfolding bandwidth is invalid, or the unfolded spacings contain degenerate or nonfinite values.
- `dCor` is `NaN` if there are fewer than three paired observations or either set of observations has zero distance variance.
