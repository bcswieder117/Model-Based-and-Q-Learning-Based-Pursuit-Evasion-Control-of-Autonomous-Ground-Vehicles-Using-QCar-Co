# Pursuit–Evasion (PEV) MATLAB Simulations

This folder contains the MATLAB scripts used to generate the numerical results, figures, CSV tables, and reproducibility manifests for Chapters 4–6 of a finite-horizon pursuit–evasion thesis study.

The simulations use an eight-state, two-vehicle kinematic-bicycle formulation:

$$
x =
\begin{bmatrix}
x_P & y_P & \psi_P & v_P & x_E & y_E & \psi_E & v_E
\end{bmatrix}^{T},
$$

where \(P\) denotes the pursuer and \(E\) denotes the evader. The controller is synthesized from a local discrete-time linearization and evaluated on nonlinear, actuator-limited kinematic-bicycle rollouts using fourth-order Runge–Kutta integration.

---

## Included Scripts

| Script | Purpose | Main Outputs |
|---|---|---|
| `PEV Ch 4 v2 with comments(1).m` | Chapter 4 model-based finite-horizon pursuit–evasion study. Computes a Riccati saddle policy and evaluates straight, limited-saddle, and full-saddle evader behaviors. | PNG figures, a nominal-results CSV summary, and a run manifest. |
| `PEV Ch 5 Code v4 with Comments(2).m` | Chapter 5 fitted-Q recovery study. Recovers stage-dependent finite-horizon saddle gains from sampled nominal transitions and compares them with the Riccati reference. | PNG diagnostics, a fitted-Q diagnostics CSV, and a run manifest. |
| `PEV Ch6 Code v3_with comments.m` | Chapter 6 robustness study. Applies wheelbase, gain, delay, actuator-limit, and noise perturbations to the fixed nominal Riccati saddle policy. | PNG figures, CSV files, and a run manifest. |

---

## Requirements

- MATLAB with standard numerical, graphics, and table-writing support.
- No external datasets or user-defined functions are required.
- The scripts are self-contained and do **not** need to be run in sequence.

The scripts use standard MATLAB functions including:

```matlab
expm
eig
rcond
cond
table
writetable
exportgraphics
