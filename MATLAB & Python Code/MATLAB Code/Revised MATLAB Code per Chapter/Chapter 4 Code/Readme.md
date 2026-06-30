
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

```matlab
expm
eig
rcond
cond
table
writetable
print
```

---

## How to Run

1. Download this folder.
2. Open MATLAB.
3. Set the MATLAB Current Folder to the directory containing the script.
4. Run the following command in the MATLAB Command Window:

```matlab
run('PEV Ch 4 v2 with comments(2).m')
```

The script automatically:

- Clears the workspace and closes existing figures.
- Initializes a fixed random seed.
- Computes the finite-horizon Riccati gains.
- Simulates all evader-policy cases.
- Creates an output directory next to the script.
- Exports figures, a CSV summary, and a run manifest.

---

## Pursuit–Evasion State Model

The two-vehicle state vector is

$$
x =
\begin{bmatrix}
x_P &
y_P &
\psi_P &
v_P &
x_E &
y_E &
\psi_E &
v_E
\end{bmatrix}^{T},
$$

where:

| Symbol | Description |
|---|---|
| $x_P, y_P$ | Pursuer planar position |
| $\psi_P$ | Pursuer heading angle |
| $v_P$ | Pursuer speed |
| $x_E, y_E$ | Evader planar position |
| $\psi_E$ | Evader heading angle |
| $v_E$ | Evader speed |

The pursuer and evader each have two control inputs:

$$
u_P =
\begin{bmatrix}
a_P \\
\delta_P
\end{bmatrix},
\qquad
u_E =
\begin{bmatrix}
a_E \\
\delta_E
\end{bmatrix},
$$

where $a$ is longitudinal acceleration and $\delta$ is steering angle.

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
| Initial pursuer position | `[-2.00, 0.00] m` |
| Initial evader position | `[0.00, 0.60] m` |
| Initial pursuer speed | `0.75 m/s` |
| Initial evader speed | `0.60 m/s` |
| Pursuer maximum acceleration | `±1.25` |
| Evader maximum acceleration | `±0.75` |
| Pursuer steering limit | `±27°` |
| Evader steering limit | `±22°` |
| Pursuer speed range | `0.00–1.45 m/s` |
| Evader speed range | `0.00–1.05 m/s` |

---

## Control Design

The controller is synthesized from a discrete-time local linearization of the two-vehicle kinematic-bicycle model.

The state penalty is defined in relative vehicle coordinates:

$$
C_{\mathrm{rel}} =
\begin{bmatrix}
I_4 & -I_4
\end{bmatrix},
$$

such that the game penalizes differences in pursuer and evader position, heading, and speed rather than tracking an arbitrary global reference point.

The nominal relative-state weighting is

$$
Q_r =
\operatorname{diag}(35,35,5,1).
$$

The base pursuer and evader control penalties are

$$
R_P =
\operatorname{diag}(0.35,0.18),
\qquad
R_E =
\operatorname{diag}(8,8).
$$

Before selecting the game parameters, the script searches through candidate state, terminal, pursuer, and evader weighting scales. A candidate is accepted only when the finite-horizon saddle-point conditions are satisfied at every time step.

These checks include:

- Positive pursuer curvature.
- Negative evader Schur-complement curvature.
- Sufficient numerical conditioning of the joint saddle matrix.

---

## Evaluated Evader Policies

The same model-based pursuer policy is evaluated against three evader behaviors.

| Evader Policy | Description |
|---|---|
| Straight | The evader applies zero control input. |
| Limited Saddle | The evader applies a scaled saddle-policy response using `evaderPolicyScale = 0.55`. |
| Full Saddle | The evader applies the complete Riccati saddle-policy response. |

This structure isolates how increasingly adversarial evader behavior affects capture time, separation, control effort, and pursuer saturation.

---

## Nonlinear Rollout

Although the controller is synthesized from a local linear model, the closed-loop simulation is evaluated using nonlinear kinematic-bicycle dynamics.

For each vehicle,

$$
\dot{x} = v\cos(\psi),
$$

$$
\dot{y} = v\sin(\psi),
$$

$$
\dot{\psi} = \frac{v}{L}\tan(\delta),
$$

$$
\dot{v} = a.
$$

The nonlinear dynamics are propagated using fourth-order Runge–Kutta integration.

Before propagation:

- Acceleration commands are clipped to vehicle limits.
- Steering commands are clipped to steering limits.
- Vehicle speeds are constrained to their permitted ranges.
- Heading angles are wrapped to the interval $[-\pi,\pi]$.

Capture occurs when the planar pursuer–evader separation satisfies

$$
\left\|
\begin{bmatrix}
x_P \\
y_P
\end{bmatrix}
-
\begin{bmatrix}
x_E \\
y_E
\end{bmatrix}
\right\|
\leq 0.35\ \mathrm{m}.
$$

---

## Output Files

The script automatically creates the following directory:

```text
PEV_Ch4_Model_Based_Output_v11/
```

The directory contains:

```text
PEV_Ch4_Model_Based_Output_v11/
├── fig4_1_trajectory_straight.png
├── fig4_2_trajectory_limited_saddle.png
├── fig4_3_trajectory_full_saddle.png
├── fig4_4_capture_time.png
├── fig4_5_capture_margin.png
├── fig4_6_pursuer_saturation.png
├── chapter4_model_based_summary.csv
└── run_manifest.txt
```

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

The file

```text
chapter4_model_based_summary.csv
```

contains the numerical data behind the Chapter 4 figures and tables.

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

The script is configured to support repeatable figure and table generation.

- The random-number generator is initialized with:

```matlab
rng(11,'twister')
```

- All simulation parameters are explicitly defined inside the script.
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
- Saturation statistics are included to identify where the unconstrained local-game commands differ from feasible vehicle commands.

---

## Thesis Citation

```text
Blaine Schwieder,
Model-Based and Q-Learning-Based Pursuit–Evasion Control of Autonomous Ground Vehicles Using QCar,
Master of Science Thesis,
Tennessee Technological University.
```
