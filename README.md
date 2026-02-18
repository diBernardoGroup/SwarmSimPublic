<img src="https://github.com/diBernardoGroup/SwarmSimPublic/assets/13195814/58bb7b7e-1966-4d4e-80ec-92e0e7fa2669" width="180">

# SwarmSim

**Simulate 2D and 3D swarms of mobile agents** with various dynamics and implement **distributed control laws** to obtain different behaviours.

To execute your first simulation simply run `Launcher.m`.

## Features
  - Simulation of swarms of mobile agents (see `Launcher.m`).
  - Implement your own dynamical model or use one of the embedded ones (see `integrateAgents.m`).
  - Implement your own distributed control law or use one of the default ones (see `globalInteractionForce.m`).
  - Acquisition of metrics to evaluate the performance.
  - Extensive simulations to study stochastic effects or different initial conditions (see `MultiLauncher.m`).
  - Extensive simulations to study the effects of the parameters (see `SequentialLauncher.m`).
  - Simulate swarms of photo-sensistive microorganisms (see `LauncherMicroorg.m`)
  - Analize [DOME](https://theopendome.org/) experiments, generate digital twins of microorganisms and run virtual experiments (see `DOME` folder) 
  - Local stability analysis of lattice configurations via linearization (see `CrystalStability.m`).
  - Interface with Robotarium code to perform advanced simulations and real life experiments (see `RobotariumSimulator.m`). 
    For more details refer to https://www.robotarium.gatech.edu.
  - Plots to visualise the simulation and the metrics.
  - Automatic save the results of the simulations.

  <img src="https://github.com/user-attachments/assets/bb57d7d6-e805-4a75-8567-b29126084fb6" width="31%">
  <img src="https://github.com/user-attachments/assets/9979db08-9035-47fa-a7d3-9d825f3fc8d3" width="30.8%">
  <img src="https://github.com/user-attachments/assets/cf3c8c66-2634-4057-b042-2f3e1995fec1" width="31%">

## Versions & Releases
**Copyrights**: If you use this code for research purposes and want to mention it in one of your publications, please cite [Giusti2023B].
  - [V3.0](https://github.com/diBernardoGroup/SwarmSimPublic/releases/tag/v3) - Authors: Andrea Giusti and Davide Salzano. Date: 2025
  - [V2.0](https://github.com/diBernardoGroup/SwarmSimPublic/releases/tag/v2) - Authors: Andrea Giusti. Date: 2023
  - [V1.1](https://github.com/diBernardoGroup/SwarmSimPublic/releases/tag/v1.1) - Authors: Andrea Giusti and Gian Carlo Maffettone. Date: 2023
  - [V1.0](https://github.com/diBernardoGroup/SwarmSimPublic/releases/tag/v1) - Authors: Andrea Giusti and Gian Carlo Maffettone. Date: 2022

## External resources
The [DOME](https://theopendome.org/) is an open-source platform for the control of microscale agents using light. 
To use the DOME, perform experiments and analyze the resulting data use [DOME-software](https://github.com/andreagiusti96/DOME-software).

The [Robotarium](https://www.robotarium.gatech.edu) allows to remotely conduct real-life swarm robotics experiments.

For additional information contact andrea.giusti@unina.it.

## References
This project implements the algorithms described in the following works:

  - \[[Giusti2026](https://royalsocietypublishing.org/rsif/article/23/234/20250780/479767)\] Giusti, A., Salzano, D., di Bernardo, M., & Gorochowski, T.E. (2026). *Data-driven inference of digital twins for high-throughput phenotyping of motile and light-responsive microorganisms*. Journal of the Royal Society Interface.
  - \[[Giusti2024](https://doi.org/10.48550/arXiv.2501.00110)\] Giusti, A. (2024). *Modelling and Control of Spatial Behaviours in Multi-Agent Systems with Applications to Biology and Robotics*. PhD Thesis.
  - \[[Giusti2023B](https://www.frontiersin.org/journals/robotics-and-ai/articles/10.3389/frobt.2023.1219931)\] Giusti, A., Maffettone, G. C., Fiore, D., Coraggio, M., & di Bernardo, M. (2023). *Distributed Control for Geometric Pattern Formation of Large-Scale Multirobot Systems*. Frontiers in Robotics and AI.
  - \[[Giusti2023A](https://ieeexplore.ieee.org/abstract/document/10160116)\] Giusti, A., Coraggio, M., & di Bernardo, M. (2023). *Local Convergence of Multi-Agent Systems Toward Rigid Lattices*. IEEE Control Systems Letters.

For [Giusti2023B] see [SwarmSimV1](https://github.com/diBernardoGroup/SwarmSimPublic/tree/SwarmSimV1), the Media folder contains the related supplementary video.
