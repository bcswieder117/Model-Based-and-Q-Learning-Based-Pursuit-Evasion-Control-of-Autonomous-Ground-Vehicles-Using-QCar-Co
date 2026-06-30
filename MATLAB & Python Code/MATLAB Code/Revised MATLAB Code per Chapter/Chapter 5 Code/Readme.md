# Chapter 5: Stage-Wise Fitted-Q Recovery for Pursuit–Evasion Control

This folder contains the MATLAB implementation for the Chapter 5 fitted-Q recovery study.

The script uses sampled transitions from the nominal finite-horizon pursuit–evasion game to recover stage-dependent quadratic Q-functions and feedback gains. The recovered gains are compared with the analytical finite-horizon Riccati saddle-policy gains used as the model-based reference.

---

## Script

| File | Description |
|---|---|
| `PEV Ch 5 Code v4 with Comments(3).m` | Recovers finite-horizon Riccati saddle gains using stage-wise fitted quadratic Q-learning and exports Chapter 5 diagnostics. |

---

## Requirements

- MATLAB
- Standard MATLAB numerical, graphics, and table-writing functionality
- No external datasets are required
- No user-defined helper functions are required

The script uses standard MATLAB functions such as:

```matlab
expm
eig
rank
cond
rcond
semilogy
table
writetable
print
```

---

## How to Run

1. Open MATLAB.
2. Set the MATLAB Current Folder to the directory containing the script.
3. Run the following command in the MATLAB Command Window:

```matlab
run('PEV Ch 5 Code v4 with Comments(3).m')
```

The script automatically:

- Clears the workspace and closes existing figures.
- Initializes a fixed random seed.
- Reconstructs the nominal finite-horizon Riccati saddle game.
- Generates sampled transitions from the nominal discrete linear model.
- Recovers fitted-Q gains at every stage.
- Evaluates Riccati and fitted-Q policies on the same nonlinear rollout.
- Creates an output directory next to the script.

---

## Study Purpose

Chapter 5 asks one focused question:

> Can a stage-wise fitted quadratic Q-function recover the finite-horizon Riccati saddle policy from sampled nominal transitions?

The Riccati solution is retained as an analytical reference. The fitted-Q method is not treated as a separate controller design problem with a changed objective, changed model, or different initial condition.

Both policies are evaluated under the same:

- Initial engagement.
- Finite horizon.
- Vehicle limits.
- Nonlinear kinematic-bicycle rollout.
- Full-saddle evader behavior.
- Sampled capture condition.

---

## Pursuit–Evasion State Model

The two-vehicle state vector is:

```text
x = [xP yP psiP vP xE yE psiE vE]'
```

where:

| Symbol | Description |
|---|---|
| `xP`, `yP` | Pursuer planar position |
| `psiP` | Pursuer heading angle |
| `vP` | Pursuer speed |
| `xE`, `yE` | Evader planar position |
| `psiE` | Evader heading angle |
| `vE` | Evader speed |

Each vehicle has longitudinal acceleration and steering-angle control inputs:

```text
uP = [aP deltaP]'
uE = [aE deltaE]'
```

The model is the same nominal two-vehicle kinematic-bicycle construction used for the Chapter 4 model-based benchmark.

---

## Simulation Configuration

| Quantity | Value |
|---|---:|
| Sampling interval | `Ts = 0.02 s` |
| Finite horizon | `N = 400` steps |
| Total horizon duration | `8.0 s` |
| Pursuer wheelbase | `Lp = 0.256 m` |
| Evader wheelbase | `Le = 0.256 m` |
| Capture radius | `0.35 m` |
| Initial pursuer state | `[-2.00, 0.00, 0.00, 0.75]` |
| Initial evader state | `[0.00, 0.60, 0.00, 0.60]` |
| Pursuer acceleration limit | `+/-1.25` |
| Evader acceleration limit | `+/-0.75` |
| Pursuer steering limit | `+/-27 degrees` |
| Evader steering limit | `+/-22 degrees` |
| Pursuer speed range | `0.00 to 1.45 m/s` |
| Evader speed range | `0.00 to 1.05 m/s` |
| Random-number generator | `rng(11,'twister')` |

---

## Nominal Riccati Reference

The script first reconstructs the nominal finite-horizon zero-sum Riccati saddle game.

The state cost is formulated in relative vehicle coordinates, so the game penalizes pursuer–evader differences in position, heading, and speed.

The nominal relative-state weighting is:

```text
Qr = diag([35, 35, 5, 1])
```

The base pursuer and evader input penalties are:

```text
Rp = diag([0.35, 0.18])
Re = diag([8.0, 8.0])
```

The script searches through candidate state, terminal, pursuer, and evader weighting scales. A candidate game is retained only when the stage-wise saddle conditions remain valid for the complete finite horizon.

The checks include:

- Positive pursuer curvature.
- Negative evader Schur-complement curvature.
- Sufficient conditioning of the joint saddle matrix.

The resulting Riccati gains serve only as the reference policy for fitted-Q recovery.

---

## Fitted-Q Recovery Method

The fitted-Q method uses a symmetric quadratic Q-function in the augmented state-action vector:

```text
zeta = [x' uP' uE']'
```

The augmented vector has:

```text
8 state variables
2 pursuer inputs
2 evader inputs
12 total variables
```

A symmetric quadratic form with 12 variables contains:

```text
78 unique quadratic coefficients
```

At each finite-horizon stage, the script:

1. Draws sampled state and action tuples.
2. Propagates each tuple through the nominal discrete linear model.
3. Builds Bellman targets using the stage cost and next-stage value matrix.
4. Constructs the quadratic feature matrix.
5. Solves a regularized least-squares problem.
6. Rebuilds the symmetric fitted Q-function matrix.
7. Checks the fitted Q-function curvature conditions.
8. Extracts fitted pursuer and evader feedback gains.
9. Propagates the recovered value matrix backward one stage.

The learning procedure proceeds backward from the terminal objective, matching the finite-horizon structure of the Riccati recursion.

---

## Sampling and Regression Settings

| Quantity | Setting |
|---|---:|
| Augmented state-action dimension | `12` |
| Unique quadratic coefficients | `78` |
| Samples per stage | `max(850, 7 × 78)` |
| Nominal samples per stage | `850` |
| Ridge regularization | `lambda = 1e-10` |

The feature matrix is checked at each stage to ensure that it has full rank before the least-squares solution is accepted.

---

## Recovery Diagnostics

The script records several diagnostics used to evaluate fitted-Q recovery.

| Diagnostic | Description |
|---|---|
| Relative Bellman residual | In-sample regression consistency of the fitted quadratic Q-function. |
| Feature rank | Confirms whether the sampled quadratic feature matrix is identifiable. |
| Feature condition | Indicates numerical sensitivity of the least-squares system. |
| Pursuer gain error | Relative difference between fitted-Q and Riccati pursuer gains. |
| Evader gain error | Relative difference between fitted-Q and Riccati evader gains. |
| Nonlinear separation history | Compares full-saddle Riccati and fitted-Q behavior on a shared nonlinear rollout. |

---

## Nonlinear Closed-Loop Recovery Check

After fitted-Q recovery, the script evaluates two controllers:

| Controller | Description |
|---|---|
| Riccati saddle policy | Exact finite-horizon model-based reference gain sequence. |
| Fitted-Q recovered policy | Gain sequence extracted from the learned stage-wise quadratic Q-functions. |

The only difference between the two rollout cases is the gain sequence. The following items remain identical:

- Initial state.
- Vehicle limits.
- Nonlinear kinematic-bicycle propagation.
- Fourth-order Runge–Kutta integration.
- Full-saddle evader response.
- Capture radius.
- Finite-horizon duration.

For each vehicle, the nonlinear rollout uses:

```text
dx/dt   = v cos(psi)
dy/dt   = v sin(psi)
dpsi/dt = (v/L) tan(delta)
dv/dt   = a
```

Acceleration and steering inputs are saturated before propagation. Vehicle speeds are constrained to their specified ranges, and heading angles are wrapped to `[-pi, pi]`.

Capture occurs when pursuer–evader planar separation is less than or equal to `0.35 m`.

---

## Output Files

The script automatically creates:

```text
PEV_Ch5_FittedQ_Recovery_Output_v11/
```

The output directory contains:

```text
PEV_Ch5_FittedQ_Recovery_Output_v11/
├── fig5_1_bellman_residuals.png
├── fig5_2_gain_recovery_error.png
├── fig5_3_full_saddle_separation.png
├── chapter5_fittedq_diagnostics.csv
└── run_manifest.txt
```

---

## Generated Figures

| Figure | Description |
|---|---|
| `fig5_1_bellman_residuals.png` | Relative Bellman residual across the finite-horizon stages. |
| `fig5_2_gain_recovery_error.png` | Relative pursuer and evader gain-recovery error between fitted-Q and Riccati gains. |
| `fig5_3_full_saddle_separation.png` | Pursuer–evader separation histories under the Riccati and fitted-Q recovered policies. |

All figures are exported as PNG files at 300 dpi.

---

## CSV Diagnostics

The file `chapter5_fittedq_diagnostics.csv` contains the numerical results behind the Chapter 5 figures and tables.

| Column | Description |
|---|---|
| `UniqueQuadraticCoefficients` | Number of independent quadratic Q-function coefficients. |
| `SamplesPerStage` | Number of sampled transitions used at each stage. |
| `MinimumFeatureRank` | Lowest feature-matrix rank observed across all stages. |
| `MaximumFeatureCondition` | Largest regularized feature-system condition value. |
| `RelativePursuerGainError` | Aggregate relative fitted-Q pursuer-gain error. |
| `RelativeEvaderGainError` | Aggregate relative fitted-Q evader-gain error. |
| `MaximumBellmanResidual` | Largest relative in-sample Bellman residual. |
| `RiccatiCaptureTime_s` | Full-saddle capture time under the Riccati policy. |
| `FittedQCaptureTime_s` | Full-saddle capture time under the fitted-Q policy. |
| `RiccatiMinimumSeparation_m` | Minimum separation under the Riccati policy. |
| `FittedQMinimumSeparation_m` | Minimum separation under the fitted-Q policy. |

---

## Reproducibility

- The script sets the random-number generator with `rng(11,'twister')`.
- All model, cost, sampling, and rollout parameters are defined inside the script.
- The fitted-Q regression uses the same reproducible sequence of sampled transitions on each execution.
- The output folder is created automatically.
- Re-running the script overwrites output files with the same names.
- The `run_manifest.txt` file records the main recovery-study assumptions.

---

## Scope of Results

This script is a matched nominal fitted-Q policy-recovery study.

The results should be interpreted within the following scope:

- Training transitions are sampled from the nominal discrete linear game.
- The fitted-Q method is evaluated against the Riccati saddle solution for the same finite-horizon objective.
- The Bellman residual is an in-sample recovery diagnostic, not a held-out generalization result.
- The nonlinear rollout is a common-policy comparison, not an independent model-based-versus-reinforcement-learning competition.
- The simulation does not represent completed physical QCar hardware validation.

---

## Thesis Citation

```text
Blaine Schwieder
Model-Based and Q-Learning-Based Pursuit-Evasion Control of Autonomous Ground Vehicles Using QCar-Compatible State-Space Models
Master of Science Thesis
Tennessee Technological University
```
