# Chapter 4: Finite-Horizon Model-Based Pursuit–Evasion Control

This folder contains the MATLAB implementation for the Chapter 4 model-based pursuit–evasion study.

The script computes a finite-horizon zero-sum Riccati saddle policy for a pursuer–evader game, then evaluates the resulting controller on nonlinear, actuator-limited kinematic-bicycle vehicle dynamics.

---

## Script

| File | Description |
|---|---|
| `PEV Ch 4 v2 with comments(2).m` | Computes the finite-horizon Riccati saddle policy, evaluates three evader behaviors, and exports Chapter 4 figures and numerical results. |

---

## Requirements

- MATLAB
- Standard MATLAB numerical and graphics functionality
- No external datasets are required
- No user-defined helper functions are required

The script uses standard MATLAB functions such as:

    expm
    eig
    rcond
    cond
    table
    writetable
    print

---

## How to Run

1. Open MATLAB.
2. Set the MATLAB Current Folder to the directory containing the script.
3. Run the following command in the MATLAB Command Window:

    run('PEV Ch 4 v2 with comments(2).m')

The script automatically clears the workspace, sets a fixed random seed, computes the finite-horizon Riccati gains, simulates all evader-policy cases, and creates an output directory next to the script.

---

## Pursuit–Evasion Model

The two-vehicle state vector is:

    x = [xP yP psiP vP xE yE psiE vE]'

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

    uP = [aP deltaP]'
    uE = [aE deltaE]'

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

---

## Control Design

The controller is synthesized from a discrete-time local linearization of the two-vehicle kinematic-bicycle model.

The cost is formulated in relative vehicle coordinates, so the game penalizes pursuer–evader differences in position, heading, and speed rather than tracking a fixed global reference point.

The nominal relative-state weighting is:

    Qr = diag([35, 35, 5, 1])

The base pursuer and evader input penalties are:

    Rp = diag([0.35, 0.18])
    Re = diag([8.0, 8.0])

Before selecting the final game weights, the script searches through candidate state, terminal, pursuer, and evader weighting scales. A candidate is accepted only when the stage-wise finite-horizon saddle conditions are satisfied.

The checks include:

- Positive pursuer curvature
- Negative evader Schur-complement curvature
- Sufficient numerical conditioning of the joint saddle matrix

---

## Evaluated Evader Policies

The same model-based pursuer policy is evaluated against three evader behaviors.

| Evader Policy | Description |
|---|---|
| Straight | The evader applies zero control input. |
| Limited Saddle | The evader applies a scaled saddle-policy response using `evaderPolicyScale = 0.55`. |
| Full Saddle | The evader applies the complete Riccati saddle-policy response. |

This setup isolates how increasingly adversarial evader behavior affects capture time, separation, control effort, and pursuer saturation.

---

## Nonlinear Rollout

Although the controller is synthesized from a local linear model, the closed-loop simulation is evaluated using nonlinear kinematic-bicycle dynamics.

For each vehicle:

    dx/dt   = v cos(psi)
    dy/dt   = v sin(psi)
    dpsi/dt = (v/L) tan(delta)
    dv/dt   = a

The nonlinear dynamics are propagated using fourth-order Runge–Kutta integration.

Before propagation:

- Acceleration commands are clipped to vehicle limits.
- Steering commands are clipped to steering limits.
- Vehicle speeds are constrained to their permitted ranges.
- Heading angles are wrapped to the interval `[-pi, pi]`.

Capture occurs when the planar pursuer–evader separation is less than or equal to `0.35 m`.

---

## Output Files

The script automatically creates:

    PEV_Ch4_Model_Based_Output_v11/

The output directory contains:

    PEV_Ch4_Model_Based_Output_v11/
    ├── fig4_1_trajectory_straight.png
    ├── fig4_2_trajectory_limited_saddle.png
    ├── fig4_3_trajectory_full_saddle.png
    ├── fig4_4_capture_time.png
    ├── fig4_5_capture_margin.png
    ├── fig4_6_pursuer_saturation.png
    ├── chapter4_model_based_summary.csv
    └── run_manifest.txt

---

## Generated Figures

| Figure | Description |
|---|---|
| `fig4_1_trajectory_straight.png` | Pursuer and straight-evader nonlinear trajectories. |
| `fig4_2_trajectory_limited_saddle.png` | Pursuer and limited-saddle evader nonlinear trajectories. |
| `fig4_3_trajectory_full_saddle.png` | Pursuer and full-saddle evader nonlinear trajectories. |
| `fig4_4_capture_time.png` | Capture-time comparison across the three evader policies. |
| `fig4_5_capture_margin.png` | Sampled capture-margin comparison across evader policies. |
| `fig4_6_pursuer_saturation.png` | Pursuer input-saturation percentage across evader policies. |

All figures are exported as PNG files at 300 dpi.

---

## CSV Summary

The file `chapter4_model_based_summary.csv` contains the numerical data behind the Chapter 4 figures and tables.

| Column | Description |
|---|---|
| `EvaderPolicy` | Straight, Limited Saddle, or Full Saddle. |
| `Captured` | Indicates whether sampled capture occurred within the horizon. |
| `CaptureTime_s` | Time of first sampled capture. |
| `MinimumSeparation_m` | Minimum planar pursuer–evader separation. |
| `SampledCaptureMargin_mm` | Capture-radius margin at minimum separation. |
| `PursuerSaturation_percent` | Percentage of pre-capture steps where pursuer input saturation occurred. |
| `PursuerEffort` | Integrated squared pursuer control effort. |
| `EvaderEffort` | Integrated squared evader control effort. |
| `GameCost` | Accumulated zero-sum game cost. |

---

## Reproducibility

- The script sets the random-number generator with `rng(11,'twister')`.
- All simulation parameters are defined inside the script.
- The output folder is created automatically.
- Re-running the script overwrites output files with the same names.
- The `run_manifest.txt` file records the principal simulation settings.

---

## Scope of Results

This script evaluates a model-based finite-horizon saddle policy under constrained nonlinear kinematic-bicycle dynamics.

The results should be interpreted as a simulation-based model-based benchmark.

- The controller is synthesized around a fixed nominal operating point.
- The nonlinear rollout includes bounded acceleration, steering, and speed.
- The results do not represent completed physical QCar hardware validation.
- Saturation statistics identify where the unconstrained local-game commands differ from feasible vehicle commands.

---

## Thesis Citation

    Blaine Schwieder
    Model-Based and Q-Learning-Based Pursuit-Evasion Control of Autonomous Ground Vehicles Using QCar-Compatible State-Space Models
    Master of Science Thesis
    Tennessee Technological University
