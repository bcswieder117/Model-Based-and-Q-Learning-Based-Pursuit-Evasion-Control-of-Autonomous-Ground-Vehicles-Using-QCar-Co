# Literature Review and Research Foundations

This directory documents the literature base for the thesis:

**Model-Based and Q-Learning-Based Pursuit–Evasion Control of Autonomous Ground Vehicles Using QCar-Compatible State-Space Models**

The thesis brings together research from autonomous driving, ground-vehicle modeling, pursuit–evasion differential games, linear-quadratic dynamic games, reinforcement learning, robust control, and QCar/QLabs implementation.

The literature is used to justify the vehicle model, formulate the pursuit–evasion game, support the fitted-Q recovery method, define robustness experiments, and establish a realistic path toward laboratory-scale autonomous vehicle implementation.

---

## Purpose of This Literature Collection

The literature review supports one central objective:

> Build a shared, physically meaningful pursuit–evasion benchmark for two autonomous ground vehicles and compare a model-based finite-horizon Riccati saddle policy with a data-driven fitted-Q policy-recovery method.

The literature is not included only as background material. Each source category supports a specific modeling, control, learning, robustness, or implementation decision used in the thesis.

---

## Literature Map

| Literature Area | Role in the Thesis |
|---|---|
| Intelligent transportation systems and autonomous driving | Motivates autonomous driving as a multi-agent decision and control problem rather than only a perception problem. |
| Ground-vehicle state-space models | Supports the selection of a kinematic bicycle model with position, heading, speed, acceleration, and steering behavior. |
| Kinematic and dynamic bicycle models | Establishes the tradeoff between model tractability and vehicle realism. |
| Pursuit–evasion differential games | Provides the adversarial interaction structure: the pursuer seeks capture while the evader seeks to avoid capture. |
| Linear-quadratic dynamic games | Supports the finite-horizon Riccati saddle-point formulation used for the model-based controller. |
| Reinforcement learning and fitted-Q methods | Supports learning from sampled state transitions without direct use of analytical system matrices. |
| Safe and robust control | Motivates perturbation studies involving mismatch, delay, saturation, and noise. |
| QCar and QLabs sources | Connects the simulation model and controller inputs to a laboratory-scale autonomous vehicle platform. |
| Graph-based multi-agent pursuit–evasion | Provides a future path from the two-vehicle benchmark toward coordinated multi-agent teams. |

---

## Literature Categories

### 1. Intelligent Transportation Systems and Autonomous Driving

This group motivates the application domain.

Autonomous driving is not only a perception problem. Vehicles must make decisions in environments containing other vehicles, pedestrians, infrastructure, uncertainty, and conflicting objectives. Lane changes, merges, following, gap acceptance, and collision avoidance all involve interactions with other agents.

This literature supports the use of pursuit–evasion as a controlled benchmark for studying competitive vehicle interaction.

Key themes include:

- Connected and automated vehicles.
- Multi-agent decision making.
- Autonomous driving control.
- Reinforcement learning for transportation systems.
- Interaction under uncertainty.
- Safety-critical vehicle behavior.

---

### 2. Ground-Vehicle State-Space Models

This group supports the vehicle-model selection process.

The thesis compares several model families before selecting the kinematic bicycle model as the primary state-space representation.

| Model Class | Typical States | Typical Inputs | Thesis Role |
|---|---|---|---|
| Point-mass or double integrator | Position and velocity | Planar acceleration | Useful for basic LQ analysis but too abstract for car-like steering behavior. |
| Unicycle | Position and heading | Forward speed and angular velocity | Useful for mobile robots but not a direct front-steering vehicle model. |
| Differential drive | Pose and wheel-related states | Left and right wheel speeds | Relevant to mobile robotics but less natural for QCar-style vehicle control. |
| Kinematic bicycle | Position, heading, speed, steering-related behavior | Acceleration and steering | Primary model because it balances realism, tractability, and implementation compatibility. |
| Dynamic bicycle | Position, heading, longitudinal velocity, lateral velocity, yaw rate | Steering, force, acceleration | Higher-fidelity extension for tire-slip, lateral dynamics, and aggressive maneuvers. |
| QCar or QLabs model | Platform-dependent states and measurements | Velocity, throttle, steering, motor commands | Implementation bridge toward virtual and physical experiments. |

The selected kinematic bicycle state is:

    z = [x y psi v]'

where:

| Variable | Description |
|---|---|
| `x` | Global longitudinal position |
| `y` | Global lateral position |
| `psi` | Vehicle heading angle |
| `v` | Longitudinal vehicle speed |

The corresponding controls are:

    u = [a delta]'

where `a` is longitudinal acceleration and `delta` is steering angle.

---

### 3. Kinematic Bicycle Model Literature

The kinematic bicycle model is the primary model used in the thesis because it preserves car-like motion without introducing the full complexity of tire forces, slip angles, lateral velocity, and yaw-rate dynamics.

The model supports:

- Planar vehicle movement.
- Heading-dependent motion.
- Non-holonomic constraints.
- Steering-limited maneuvers.
- Longitudinal acceleration.
- Linearization and discretization.
- Finite-horizon LQ game design.
- MATLAB simulation.
- QCar and QLabs compatibility.

The dynamic bicycle model remains important as a higher-fidelity future extension when lateral velocity, yaw rate, tire forces, and inertial effects become central to the research question.

---

### 4. Pursuit–Evasion and Differential Game Theory

This group provides the mathematical foundation for the adversarial interaction.

In a pursuit–evasion game:

- The pursuer attempts to reduce separation and enter a capture region.
- The evader attempts to avoid capture or delay capture.
- Both vehicles must obey their own dynamics and input constraints.

The thesis uses a two-player zero-sum formulation because it provides a controlled and interpretable setting for comparing model-based and data-driven control approaches.

Important concepts include:

- Capture sets.
- Relative geometry.
- Pursuer and evader policies.
- Saddle-point equilibria.
- Finite-horizon dynamic games.
- Differential games.
- Discrete-time pursuit–evasion formulations.

The two-vehicle case is treated as the basic unit of adversarial vehicle interaction before extending toward larger coordinated teams.

---

### 5. Linear-Quadratic Dynamic Games

This literature supports the Chapter 4 model-based controller.

The nonlinear kinematic bicycle model is locally linearized and discretized to form a finite-horizon game:

    X(k+1) = A X(k) + Bp uP(k) + Be uE(k)

The model-based controller is obtained through a Riccati-type saddle-point recursion.

This literature supports:

- Stage-dependent feedback gains.
- Finite-horizon pursuit–evasion control.
- Saddle-point conditions.
- Pursuer curvature conditions.
- Evader Schur-complement conditions.
- Numerical conditioning checks.
- Terminal penalties.
- Relative-state cost construction.

The Riccati policy serves as the analytical reference for the fitted-Q recovery study.

---

### 6. Reinforcement Learning and Fitted-Q Recovery

This group supports the Chapter 5 data-driven policy-recovery method.

The fitted-Q method uses sampled transitions instead of direct access to the analytical transition matrices used by the Riccati recursion.

The method learns from:

- State samples.
- Pursuer action samples.
- Evader action samples.
- Stage costs.
- Successor states.

The fitted-Q procedure reconstructs a stage-wise quadratic action-value function and extracts the pursuer and evader feedback gains.

The purpose is not to claim that fitted-Q universally outperforms the Riccati solution. Instead, the objective is to determine whether sampled nominal transitions can recover the same finite-horizon saddle policy.

Key concepts include:

- Q-learning.
- Bellman regression.
- Quadratic Q-functions.
- Least-squares fitting.
- Function approximation.
- Feature rank.
- Numerical conditioning.
- Gain-recovery error.
- Bellman residuals.

---

### 7. Safe and Robust Control Literature

This group motivates the Chapter 6 robustness studies.

A controller that is optimal for a nominal linearized model may degrade when the actual execution conditions differ from the design assumptions.

The thesis therefore evaluates the nominal Riccati saddle policy under execution-side perturbations.

| Perturbation | Physical Interpretation |
|---|---|
| Wheelbase mismatch | Geometry or model-parameter uncertainty. |
| Steering-gain mismatch | Differences between requested and realized steering response. |
| Longitudinal-gain mismatch | Differences between requested and realized acceleration response. |
| Command delay | Computation, communication, sensing, or actuator delay. |
| Actuator-limit severity | Reduced available acceleration or steering authority. |
| Process noise | Unmodeled disturbances affecting vehicle propagation. |
| Measurement noise | Imperfect state estimates supplied to the controller. |
| Combined perturbation case | Representative execution-side deviation from the nominal model. |

These studies are sensitivity analyses. They do not claim that the selected perturbation magnitudes are identified QCar hardware parameters.

---

### 8. QCar and QLabs Implementation Literature

This group connects the mathematical model to a laboratory-scale autonomous vehicle platform.

QCar and QLabs references are used to support:

- Car-like vehicle geometry.
- Steering and velocity command interfaces.
- Actuator constraints.
- Virtual simulation workflows.
- MATLAB, Simulink, and Python implementation pathways.
- The distinction between simulation validation and completed hardware validation.

The thesis does not claim completed physical QCar validation. Instead, it develops a QCar-oriented simulation and robustness pipeline intended to reduce the gap between analytical control design and future platform testing.

---

### 9. Graph-Based Multi-Agent Pursuit–Evasion

This group supports future extensions beyond one pursuer and one evader.

The present thesis focuses on a two-vehicle interaction because it provides:

- A clear state-space model.
- A transparent game formulation.
- An analytical model-based baseline.
- A controlled fitted-Q recovery problem.
- Interpretable robustness experiments.

Graph-based multi-agent pursuit–evasion provides a natural future direction in which several pursuers and evaders interact through communication, neighbor relationships, coordination policies, and network structure.

---

## Literature-to-Thesis Traceability

| Thesis Component | Supporting Literature | Purpose |
|---|---|---|
| Chapter 1: Motivation and problem framing | Autonomous driving, ITS, multi-agent control | Establishes why vehicle interaction is a meaningful control problem. |
| Chapter 2: Literature review | All categories in this directory | Establishes the research gap and the literature-supported research direction. |
| Chapter 3: Vehicle modeling | Kinematic bicycle, dynamic bicycle, QCar, QLabs | Justifies the selected two-vehicle state-space model. |
| Chapter 4: Model-based controller | Pursuit–evasion, LQ games, Riccati theory | Supports the finite-horizon saddle-point controller. |
| Chapter 5: Fitted-Q recovery | Q-learning, Bellman regression, least-squares methods | Supports data-driven recovery from sampled transitions. |
| Chapter 6: Robustness studies | Safe RL, robust control, uncertainty analysis | Supports mismatch, delay, noise, and saturation experiments. |
| Future work | Graph theory, multi-agent control, dynamic vehicle models | Supports extensions to multi-vehicle teams and higher-fidelity vehicle dynamics. |

---

## Central Research Gap

The literature identifies a gap between three research directions:

1. Pursuit–evasion studies often use abstract agent dynamics, guidance laws, or low-order motion models.
2. Autonomous-vehicle studies often use realistic vehicle models but focus on tracking, platooning, lane keeping, cooperative driving, or collision avoidance rather than adversarial pursuit–evasion.
3. Reinforcement-learning pursuit–evasion studies often use complex simulation environments but do not compare learned policies with an exact linear-quadratic game solution on the same model.

This thesis addresses that gap by developing a shared, model-centered pursuit–evasion benchmark with:

- A literature-supported kinematic bicycle vehicle model.
- Two autonomous ground vehicles.
- A finite-horizon Riccati saddle-point controller.
- A fitted-Q recovery method using sampled nominal transitions.
- Common state definitions, capture conditions, dynamics, and performance metrics.
- Nonlinear actuator-limited rollout evaluation.
- QCar-relevant perturbation studies.
- A defined pathway toward QLabs and future QCar implementation.

---

## Recommended Literature Folder Structure

    Literature/
    ├── 01_Autonomous_Driving_and_ITS/
    ├── 02_Ground_Vehicle_Models/
    ├── 03_Kinematic_and_Dynamic_Bicycle_Models/
    ├── 04_Pursuit_Evasion_and_Differential_Games/
    ├── 05_Linear_Quadratic_Dynamic_Games/
    ├── 06_Reinforcement_Learning_and_Fitted_Q/
    ├── 07_Safe_and_Robust_Control/
    ├── 08_QCar_and_QLabs/
    ├── 09_Graph_Based_Multi_Agent_Control/
    ├── references.bib
    └── literature_notes.md

---

## Citation and Source Management

All claims in the thesis should be traceable to a peer-reviewed paper, textbook, technical report, or official platform source.

Recommended practice:

- Use peer-reviewed papers for scientific claims.
- Use textbooks for foundational control and game-theory results.
- Use official Quanser and MathWorks documentation for platform and software details.
- Verify author names, publication year, venue, page numbers, and DOI information before final submission.
- Keep the BibTeX database synchronized with the PDF library.
- Record short notes explaining how each source supports a thesis decision.

Platform documentation is useful for QCar and QLabs implementation details, but it should be distinguished from peer-reviewed research evidence.

---

## Scope Notes

This literature collection supports a controlled engineering benchmark.

It does not claim to solve all autonomous-driving interaction problems, reproduce a full traffic simulator, or provide completed physical QCar validation.

The thesis focuses on:

- One pursuer and one evader.
- Reduced-order car-like vehicle dynamics.
- A fixed local linearization for the finite-horizon game.
- A fitted-Q recovery study under matched nominal transitions.
- Execution-side perturbation studies.
- A future implementation path toward QLabs and QCar hardware.

---

## Associated Thesis

    Blaine Christopher Swieder
    Model-Based and Q-Learning-Based Pursuit–Evasion Control of Autonomous Ground Vehicles Using QCar-Compatible State-Space Models
    Master of Science Thesis
    Tennessee Technological University
    August 2026
