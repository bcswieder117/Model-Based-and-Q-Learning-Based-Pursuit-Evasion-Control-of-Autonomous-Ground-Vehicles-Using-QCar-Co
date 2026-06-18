# Model-Based and Q-Learning-Based Pursuit-Evasion Control of Autonomous Ground Vehicles Using QCar-Compatible State-Space Models

A portfolio of code, experiments, and documentation for my MSc Thesis at Tennessee Tech University. My research sits at the intersection of optimal control, game-theoretic multi-agent systems, and reinforcement learning, with applications to pursuit–evasion and autonomous vehicles.

This will be presented at the Tennessee Tech Research Day 2026, the link about the event is given here: https://www.tntech.edu/research/research-day.php

The thesis defense will be the third week of July. All materials will be uploaded then and a link to ProQuest will also be provided no later than 7/31/2026

## Abstract

Autonomous ground vehicles usually operate in environments where safety-critical behavior relies on interaction with other moving agents. This thesis aims to investigate pursuit-evasion as a controlled benchmark to compare model-based game control and data-driven policy recovery for two car-like autonomous ground vehicles. A QCar-compatible kinematic bicycle model is used such that each vehicle retains physically meaningful states, which includes planar position, heading angle, and longitudinal speed, with acceleration and steering as control inputs. 

The model-based portion formulates a finite-horizon, two-player, zero-sum linear-quadratic pursuit-evasion game. The pursuer seeks to drive the relative vehicle configuration into a prescribed capture set, while the evader seeks to delay capture. A local discrete-time game is obtained by linearizing and discretizing the nonlinear kinematic bicycle dynamics. A Riccati-type saddle-point recursion the produces stage-dependent feedback policies for both players, which are then evaluated on the nonlinear actuator-limited vehicle roll-out. 


The data-driven policy recovery uses fitted quadratic Q-learning to recover the same finite-horizon saddle policy from the sampled states, actions, costs, and successor states without providing the learner with the analytical transition matrices used by the Riccati recursion. As a result, this creates a direct comparison between an analytical controller that is derived from known dynamics and a model-free recovery procedure that is operating on the same game. 

Finally, QCar-relevant perturbations are introduced, such as wheelbase mismatch, steering, and longitudinal-command mismatch, noise, delay, and actuator-limit severity. These studies provide the transfer-readiness evidence by showing how the nominal saddle policy responds to platform relevant derivations while distinguishing simulation validation from completed hardware validation. This result is a reproducible pipeline for QCar-oriented pursuit-evasion research and hardware experiments. 
