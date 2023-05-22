# SwarmSim

Simulate 2D and 3D swarms of mobile agents with various dynamics and implement distributed control laws to obtain different behaviours.

To execute your first simulation simply run Launcher.m.

Features:
  - Simulation of swarms of mobile agents.
  - Implement your own dynamical model or use one of the embedded ones (see integrateAgents.m).
  - Implement your own distributed control law or use one of the embedded ones (see globalInteractionForce.m).
  - Acquisition of metrics to evaluate the performance.
  - Extensive simulations to study stochastic effects or different initial conditions (see MultiLauncher.m).
  - Extensive simulations to study the effects of the parameters (see SequentialLauncher.m).
  - Local stability analysis through linearization (see CrystalStability).
  - Interface with Robotarium code to perform advanced simulations and real life experiments (see RobotariumSimulator). 
    For more details and to use the Robotarium refer to https://www.robotarium.gatech.edu.
  - Plots to visualise the simulation and the metrics.
  - Automatic save the results of the simulations.

This project implements the algorithms described in [Giusti2022] (see Release v1) and [Giusti2023]. In the Media folder, there is a video supplement for [Giusti2022].

For information contact andrea.giusti@unina.it.

Copyrights: If you use this code for research purposes and want to mention it in one of your publications, please cite [Giusti2022].

Versions:
  - V2 - Authors: Andrea Giusti. Date: 2023
  - V1 - Authors: Andrea Giusti and Gian Carlo Maffettone. Date: 2022

References:
  - [Giusti2023] Giusti, A., Coraggio, M., & di Bernardo, M. (2023). Local convergence of multi-agent systems towards triangular patterns. arXiv preprint arXiv:2303.11865.
  - [Giusti2022] Giusti, A., Maffettone, G. C., Fiore, D., Coraggio, M., & di Bernardo, M. (2022). Distributed control for geometric pattern formation of large-scale multirobot systems. arXiv preprint arXiv:2207.14567.
