# Chapter 6: Model-Based Robustness and Perturbation Study

This folder contains the MATLAB implementation for the Chapter 6 robustness study of the nominal finite-horizon Riccati saddle policy.

The script evaluates how a fixed model-based pursuit–evasion controller responds to execution-side deviations relevant to future QLabs or QCar implementation, including model mismatch, command delay, actuator limits, process noise, and measurement noise.

The Chapter 6 script evaluates the nominal Riccati saddle policy only. Fitted-Q recovery is intentionally isolated in the Chapter 5 script.

---

## Script

| File | Description |
|---|---|
| `PEV Ch6 Code v3_with comments(1).m` | Evaluates wheelbase mismatch, gain mismatch, delay, actuator limits, noise, and a representative combined-perturbation case using the fixed nominal Riccati saddle policy. |

---

## Requirements

- MATLAB
- Standard MATLAB numerical, graphics, and table-writing functionality
- No external datasets are required
- No user-defined helper functions are required

The script uses standard MATLAB functions such as:

    expm
    eig
    rcond
    cond
    rng
    errorbar
    table
    writetable
    print

---

## How to Run

1. Open MATLAB.
2. Set the MATLAB Current Folder to the directory containing the script.
3. Run the following command in the MATLAB Command Window:

    run('PEV Ch6 Code v3_with comments(1).m')

The script automatically:

- Clears the workspace and closes existing figures.
- Initializes repeatable random-number-generator seeds.
- Reconstructs the nominal finite-horizon Riccati saddle game.
- Runs all deterministic perturbation sweeps.
- Runs the noise Monte Carlo study.
- Runs a representative combined-perturbation trajectory.
- Creates an output directory next to the script.
- Exports figures, CSV files, and a run manifest.

---

## Study Purpose

Chapter 6 asks one focused question:

> How does a fixed nominal model-based saddle policy respond to mismatch, delay, noise, and actuator limits?

The controller is synthesized once from the nominal finite-horizon game. It is not retuned separately for each perturbation case.

This allows the study to isolate execution-side sensitivity rather than mixing robustness effects with controller redesign.

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

The nominal controller is synthesized from a local discrete-time linearization of the two-vehicle kinematic-bicycle model.

---

## Nominal Simulation Configuration

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
| Base random seed | `rng(11,'twister')` |

---

## Nominal Riccati Saddle Policy

The script reconstructs the finite-horizon zero-sum Riccati saddle policy used as the model-based benchmark.

The state cost is formulated in relative vehicle coordinates, so the game penalizes differences in position, heading, and speed between the pursuer and evader.

The nominal relative-state weighting is:

    Qr = diag([35, 35, 5, 1])

The base pursuer and evader control penalties are:

    Rp = diag([0.35, 0.18])
    Re = diag([8.0, 8.0])

The script searches through candidate state, terminal, pursuer, and evader weighting scales. A candidate is accepted only when the finite-horizon saddle conditions remain valid at every stage.

The checks include:

- Positive pursuer curvature.
- Negative evader Schur-complement curvature.
- Sufficient conditioning of the joint saddle matrix.

The selected Riccati gain sequence remains fixed for all Chapter 6 perturbation cases.

---

## Nonlinear Vehicle Rollout

Although the controller is synthesized from a local linear model, all perturbation cases are evaluated using nonlinear kinematic-bicycle propagation.

For each vehicle:

    dx/dt   = v cos(psi)
    dy/dt   = v sin(psi)
    dpsi/dt = (v/L) tan(delta)
    dv/dt   = a

The nonlinear dynamics are propagated using fourth-order Runge–Kutta integration.

Before propagation:

- Acceleration commands are clipped to actuator limits.
- Steering commands are clipped to steering limits.
- Vehicle speeds are constrained to their allowed ranges.
- Heading angles are wrapped to `[-pi, pi]`.

Capture occurs when the planar pursuer–evader separation is less than or equal to `0.35 m`.

---

## Perturbation Studies

Each deterministic sweep changes one execution-side mechanism while maintaining the same nominal Riccati gain sequence.

| Study | Tested Values | Description |
|---|---|---|
| Wheelbase mismatch | `0.90, 0.95, 1.00, 1.05, 1.10` | Scales both vehicle wheelbases during nonlinear propagation. |
| Steering-gain mismatch | `0.85, 0.95, 1.00, 1.05, 1.15` | Scales steering commands before saturation. |
| Longitudinal-gain mismatch | `0.85, 0.95, 1.00, 1.05, 1.15` | Scales acceleration commands before saturation. |
| Command delay | `0, 1, 2, 3, 4` steps | Applies delayed control commands to the nonlinear rollout. |
| Actuator-limit severity | `0.70, 0.80, 0.90, 1.00` | Scales available acceleration and steering authority. |
| Noise sensitivity | `0.000, 0.005, 0.015, 0.030` | Adds process noise and measurement noise across repeated trials. |

For the deterministic sweeps, the script records:

- Capture status.
- Capture time.
- Minimum separation.
- Pursuer acceleration saturation percentage.
- Pursuer steering saturation percentage.

---

## Noise Monte Carlo Study

The noise study uses:

| Quantity | Setting |
|---|---:|
| Noise scales | `0.000, 0.005, 0.015, 0.030` |
| Trials per noise scale | `200` |
| Bootstrap replications | `10,000` |
| Measurement-noise scale | `0.5 × process-noise scale` |
| Capture-rate uncertainty | 95% Wilson confidence interval |
| Time uncertainty | 95% percentile bootstrap confidence interval |

The script reports:

- Capture rate.
- Lower and upper 95% Wilson confidence bounds for capture rate.
- Horizon-censored mean capture time.
- Lower and upper 95% bootstrap confidence bounds for mean capture time.

For the horizon-censored mean-time calculation, a trial without capture is assigned the full simulation horizon:

    N * Ts = 8.0 s

This makes the reported time metric account for both delayed capture and failed capture within the finite horizon.

---

## Representative Combined-Perturbation Case

The script also generates one representative trajectory that combines moderate deviations:

| Quantity | Value |
|---|---:|
| Wheelbase scale | `1.05` |
| Steering-gain scale | `0.95` |
| Longitudinal-gain scale | `0.95` |
| Actuator-limit scale | `0.90` |
| Command delay | `1` step |
| Process-noise scale | `0.002` |
| Measurement-noise scale | `0.001` |
| Random seed | `rng(777,'twister')` |

This trajectory is intended as an execution-side sensitivity example for a future QLabs or QCar transition. It is not an identified hardware parameter set or a completed physical vehicle validation.

---

## Output Files

The script automatically creates:

    PEV_Ch6_Perturbation_Model_Based_Output_v10/

The output directory contains:

    PEV_Ch6_Perturbation_Model_Based_Output_v10/
    ├── fig6_1_wheelbase_capture_time_THESIS.png
    ├── fig6_2_steering_capture_time_THESIS.png
    ├── fig6_3_longitudinal_capture_time_THESIS.png
    ├── fig6_4_noise_capture_rate_THESIS.png
    ├── fig6_4b_noise_horizon_censored_time_THESIS.png
    ├── fig6_5_delay_capture_time_THESIS.png
    ├── fig6_6_actuator_limit_capture_time_THESIS.png
    ├── fig6_7_steering_saturation_THESIS.png
    ├── fig6_8_combined_trajectory_THESIS.png
    ├── fig6_9_combined_separation_THESIS.png
    ├── chapter6_combined_trajectory_THESIS.csv
    ├── chapter6_model_based_perturbation_summary_v10.csv
    ├── chapter6_model_based_noise_monte_carlo_v10.csv
    └── run_manifest.txt

All figures are exported as PNG files at 300 dpi.

---

## Generated Figures

| Figure | Description |
|---|---|
| `fig6_1_wheelbase_capture_time_THESIS.png` | Capture time versus wheelbase scale. |
| `fig6_2_steering_capture_time_THESIS.png` | Capture time versus steering-gain scale. |
| `fig6_3_longitudinal_capture_time_THESIS.png` | Capture time versus longitudinal-gain scale. |
| `fig6_4_noise_capture_rate_THESIS.png` | Capture rate under noise with 95% Wilson confidence intervals. |
| `fig6_4b_noise_horizon_censored_time_THESIS.png` | Horizon-censored mean time under noise with 95% bootstrap confidence intervals. |
| `fig6_5_delay_capture_time_THESIS.png` | Capture time versus command delay. |
| `fig6_6_actuator_limit_capture_time_THESIS.png` | Capture time versus available actuator authority. |
| `fig6_7_steering_saturation_THESIS.png` | Pursuer steering saturation versus actuator-limit scale. |
| `fig6_8_combined_trajectory_THESIS.png` | Pursuer and evader trajectories for the combined-perturbation case. |
| `fig6_9_combined_separation_THESIS.png` | Pursuer–evader separation history for the combined-perturbation case. |

---

## CSV Files

### Perturbation Summary

The file:

    chapter6_model_based_perturbation_summary_v10.csv

contains the deterministic-sweep and noise-summary results.

| Column | Description |
|---|---|
| `perturbation` | Type of perturbation study. |
| `level` | Tested perturbation value. |
| `controller` | Evaluated policy. |
| `captured` | Indicates whether capture occurred. |
| `capture_time_s` | Capture time or horizon-censored time. |
| `min_separation_m` | Minimum planar separation. |
| `pursuer_accel_saturation_percent` | Pursuer acceleration saturation percentage. |
| `pursuer_steering_saturation_percent` | Pursuer steering saturation percentage. |

### Noise Monte Carlo Results

The file:

    chapter6_model_based_noise_monte_carlo_v10.csv

contains the noise-study statistics.

| Column | Description |
|---|---|
| `noise_scale` | Process-noise scale. |
| `trials_per_condition` | Number of Monte Carlo trials. |
| `capture_rate_percent` | Percentage of trials resulting in capture. |
| `capture_rate_ci95_low` | Lower 95% Wilson confidence bound. |
| `capture_rate_ci95_high` | Upper 95% Wilson confidence bound. |
| `horizon_censored_mean_time_s` | Mean time with failed captures assigned the full horizon. |
| `mean_time_ci95_low` | Lower 95% bootstrap confidence bound. |
| `mean_time_ci95_high` | Upper 95% bootstrap confidence bound. |

### Combined-Trajectory Data

The file:

    chapter6_combined_trajectory_THESIS.csv

contains the stored trajectory for the representative combined-perturbation case.

| Column | Description |
|---|---|
| `time_s` | Simulation time. |
| `xP`, `yP`, `psiP`, `vP` | Pursuer position, heading, and speed. |
| `xE`, `yE`, `psiE`, `vE` | Evader position, heading, and speed. |
| `separation_m` | Planar pursuer–evader separation. |

---

## Reproducibility

- The script initializes the nominal run with `rng(11,'twister')`.
- Separate deterministic seeds are used for individual sweep conditions.
- Separate deterministic seeds are used for Monte Carlo trials.
- The representative combined case uses `rng(777,'twister')`.
- All model, controller, perturbation, and uncertainty settings are defined inside the script.
- The output folder is created automatically.
- Re-running the script overwrites files with the same names.
- The `run_manifest.txt` file records the main simulation, controller, and uncertainty settings.

---

## Scope of Results

This script is a model-based robustness and sensitivity study.

The results should be interpreted within the following scope:

- One nominal finite-horizon Riccati saddle policy is evaluated.
- The policy is not retrained or retuned for each perturbation condition.
- The perturbation scales are transparent sensitivity settings rather than identified QCar hardware parameters.
- The noise settings are not a calibrated QCar covariance model.
- The study is not a completed physical QCar validation.
- Fitted-Q recovery is intentionally excluded from this script and is reported separately in Chapter 5.

---

## Thesis Citation

    Blaine Schwieder
    Model-Based and Q-Learning-Based Pursuit-Evasion Control of Autonomous Ground Vehicles Using QCar-Compatible State-Space Models
    Master of Science Thesis
    Tennessee Technological University
