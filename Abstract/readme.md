# Model-Based and Q-Learning-Based Pursuit–Evasion Control of Autonomous Ground Vehicles Using QCar-Compatible State-Space Models

> **Project Status:** Active MSc thesis project at Tennessee Technological University. This repository contains the reproducible MATLAB simulations, code, literature documentation, and supporting materials developed for the thesis.

## Overview

This research develops a controlled pursuit–evasion benchmark for two autonomous ground vehicles. The project connects optimal control, differential games, data-driven policy recovery, and autonomous-vehicle modeling through a shared simulation environment.

The benchmark compares a finite-horizon model-based Riccati saddle policy with a stage-wise fitted-Q recovery method. Both approaches use the same vehicle-state definition, capture condition, initial engagement geometry, nonlinear rollout, actuator limits, and performance metrics.

## Abstract

Autonomous ground vehicles operate in environments where safe behavior depends on interaction with other moving agents. This thesis develops a controlled pursuit–evasion benchmark for comparing model-based game control with data-driven policy recovery for two car-like autonomous ground vehicles.

Each vehicle is represented by a QCar-compatible kinematic bicycle model with physically meaningful states: planar position, heading angle, and longitudinal speed. Acceleration and steering angle serve as the control inputs. This representation preserves car-like motion and steering geometry while remaining suitable for linearization, finite-horizon game design, sampled-data learning, MATLAB simulation, and future QCar-oriented implementation.

The model-based component formulates a finite-horizon, two-player, zero-sum linear-quadratic pursuit–evasion game. The pursuer seeks to drive the relative vehicle configuration into a prescribed capture set, while the evader seeks to delay or avoid capture. A local discrete-time game model is obtained by linearizing and discretizing the nonlinear kinematic bicycle dynamics. A Riccati-type saddle-point recursion then produces stage-dependent feedback policies for both players, which are evaluated on nonlinear rollouts with acceleration, steering, and speed constraints.

The data-driven component uses stage-wise fitted quadratic Q-learning to recover the same finite-horizon saddle policy from sampled states, actions, costs, and successor states without providing the learner with the analytical transition matrices used by the Riccati recursion. This creates a direct comparison between an analytical controller derived from known dynamics and a data-driven recovery procedure operating on the same pursuit–evasion game.

Finally, the nominal model-based saddle policy is evaluated under QCar-relevant perturbations, including wheelbase mismatch, steering and longitudinal-command mismatch, command delay, process noise, measurement noise, and actuator-limit severity. These experiments provide simulation-based transfer-readiness evidence while clearly distinguishing robustness evaluation from completed physical QCar hardware validation.

## Research Focus

- QCar-compatible kinematic bicycle modeling for autonomous ground vehicles
- Finite-horizon zero-sum linear-quadratic pursuit–evasion games
- Riccati saddle-point control
- Stage-wise fitted quadratic Q-learning
- Nonlinear actuator-limited vehicle rollouts
- Robustness to mismatch, delay, noise, and actuator limitations
- Future transfer toward QLabs and QCar implementation

## Scope

This repository presents a reproducible simulation and analysis pipeline. It does not claim completed physical QCar validation, universal reinforcement-learning performance, or a full autonomous-driving traffic simulator.

The current work is limited to:

- One pursuer and one evader
- Reduced-order kinematic bicycle dynamics
- A fixed nominal operating-point linearization
- A finite-horizon model-based game solution
- Matched nominal fitted-Q policy recovery
- Nonlinear actuator-limited rollout evaluation
- QCar-relevant perturbation and sensitivity studies

## Repository Contents

- MATLAB implementations for Chapters 4–6
- Model-based finite-horizon pursuit–evasion simulations
- Fitted-Q recovery diagnostics
- Robustness and perturbation experiments
- Literature documentation and thesis-supporting notes
- Reproducible figures, CSV summaries, and run manifests

## Citation

Citation information for the final archived thesis will be added after thesis submission and approval.
