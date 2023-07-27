# SwarmSimV1

This project implements the algorithms described in [Giusti2022].

Simulate swarms of planar mobile agents with first order dynamics and implement distributed control laws to obtain geometric lattice formation.

To execute your first simulation simply run Launcher.m.

Features:
  - Simulation of swarms of planar mobile agents with first order dynamics.
  - Distributed control laws to obtain geometric lattice formation.
  - Acquisition of metrics to evaluate the performance.
  - Adaptive tuning of the control gains.
  - Embedded options to test agents failure, dynamic lattice reconfiguration and actuation noise.
  - Extensive simulations to tune the control gains (see BruteForceTuning.m).
  - Extensive simulations to study the effects of the parameters (see StabilityAnalysis and SequentialLauncher.m).
  - Local stability analysis through linearization (see crystalStabilityMulti).
  - Interface with Robotarium code to perform advanced simulations and real life experiments (see RobotariumSimulator). 
    For more details and to use the Robotarium refer to https://www.robotarium.gatech.edu.
  - Plot functions to visualise the results.

For more information contact andrea.giusti@unina.it.

Copyrights: If you use this code for research purposes and want to mention it in one of your publications, please cite [Giusti2022].

Versions:
  - [V1.1](https://github.com/diBernardoGroup/SwarmSimPublic/releases/tag/v1.1) - Authors: Andrea Giusti and Gian Carlo Maffettone. Date: 2023
  - [V1.0](https://github.com/diBernardoGroup/SwarmSimPublic/releases/tag/v1) - Authors: Andrea Giusti and Gian Carlo Maffettone. Date: 2022
    
For more advanced features checkout the current version of [SwarmSim](https://github.com/diBernardoGroup/SwarmSimPublic/tree/main).

In the Media folder, there is a video supplement for [Giusti2022].

\[[Giusti2022](https://arxiv.org/abs/2207.14567)\] Giusti, A., Maffettone, G. C., Fiore, D., Coraggio, M., & di Bernardo, M. (2022). Distributed Control for Geometric Pattern Formation of Large-Scale Multirobot Systems. arXiv preprint arXiv:2207.14567.
